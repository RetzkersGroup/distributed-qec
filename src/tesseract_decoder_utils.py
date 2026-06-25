"""Sinter wrapper around the Tesseract decoder.

Provides a `Sampler` (not just a `Decoder`) so that per-shot
`low_confidence_shots` get written into the `custom_counts` column of the CSV
(the `Decoder` interface only reports shots/errors/discards). Two presets,
`tesseract_short` and `tesseract_long`, mirror the ones in
`stairway-codes/src/decoder_utils.py`.
"""

from __future__ import annotations

import collections
import time
from typing import Dict

import numpy as np
import stim
from sinter import CompiledDecoder, Decoder
from sinter._data import AnonTaskStats, Task
from sinter._decoding._stim_then_decode_sampler import (
    StimThenDecodeSampler,
    _CompiledStimThenDecodeSampler,
    classify_discards_and_errors,
)


TESSERACT_LONG_KWARGS = {
    "pqlimit": 1_000_000,
    "det_beam": 20,
    "beam_climbing": True,
    "no_revisit_dets": True,
}

TESSERACT_SHORT_KWARGS = {
    "pqlimit": 200_000,
    "det_beam": 15,
    "beam_climbing": True,
    "no_revisit_dets": True,
}


class TesseractCompiledSinterDecoder(CompiledDecoder):
    def __init__(self, decoder):
        self.decoder = decoder
        self.low_confidence_shots = 0
        self.last_batch_low_confidence_shots = 0

    def decode_shots_bit_packed(
        self, *, bit_packed_detection_event_data: np.ndarray,
    ) -> np.ndarray:
        num_shots = bit_packed_detection_event_data.shape[0]
        num_detectors = self.decoder.config.dem.num_detectors
        num_observables = self.decoder.config.dem.num_observables

        unpacked_shots = np.unpackbits(
            bit_packed_detection_event_data,
            axis=1,
            count=num_detectors,
            bitorder="little",
        ).astype(bool)

        predictions = np.zeros((num_shots, num_observables), dtype=bool)
        self.last_batch_low_confidence_shots = 0
        for i in range(num_shots):
            predictions[i] = self.decoder.decode(unpacked_shots[i])
            if self.decoder.low_confidence_flag:
                self.last_batch_low_confidence_shots += 1

        if self.last_batch_low_confidence_shots:
            print(
                f"WARNING: tesseract low confidence "
                f"{self.last_batch_low_confidence_shots}/{num_shots}"
            )

        self.low_confidence_shots += self.last_batch_low_confidence_shots
        return np.packbits(predictions, axis=1, bitorder="little")


class TesseractSinterDecoder(Decoder):
    def __init__(self, **tesseract_kwargs):
        self.tesseract_kwargs = tesseract_kwargs

    def compile_decoder_for_dem(
        self, *, dem: stim.DetectorErrorModel,
    ) -> CompiledDecoder:
        import tesseract_decoder.tesseract as tesseract

        config = tesseract.TesseractConfig(dem=dem, **self.tesseract_kwargs)
        decoder = tesseract.TesseractDecoder(config)
        return TesseractCompiledSinterDecoder(decoder)


class TesseractSampler(StimThenDecodeSampler):
    """Stim sampler that drives `TesseractSinterDecoder` and records
    `low_confidence_shots` in `custom_counts`."""

    def compiled_sampler_for_task(
        self, task: Task,
    ) -> "_CompiledTesseractSampler":
        assert isinstance(self.decoder, TesseractSinterDecoder)
        return _CompiledTesseractSampler(
            decoder=self.decoder,
            task=task,
            count_detection_events=self.count_detection_events,
            count_observable_error_combos=self.count_observable_error_combos,
            tmp_dir=self.tmp_dir,
        )


class _CompiledTesseractSampler(_CompiledStimThenDecodeSampler):
    def sample(self, max_shots: int) -> AnonTaskStats:
        t0 = time.monotonic()
        dets, actual_obs = self.stim_sampler.sample(
            shots=max_shots, bit_packed=True, separate_observables=True,
        )
        num_shots = dets.shape[0]

        custom_counts = collections.Counter()
        if self.count_detection_events:
            custom_counts["detectors_checked"] += self.num_det * num_shots
            for b in range(8):
                custom_counts["detection_events"] += np.count_nonzero(dets & (1 << b))

        if self.task.postselection_mask is not None:
            discarded_flags = np.any(dets & self.task.postselection_mask, axis=1)
            num_discards_1 = int(np.count_nonzero(discarded_flags))
            if num_discards_1:
                dets = dets[~discarded_flags, :]
                actual_obs = actual_obs[~discarded_flags, :]
        else:
            num_discards_1 = 0

        predictions = self.compiled_decoder.decode_shots_bit_packed(
            bit_packed_detection_event_data=dets,
        )

        num_discards_2, num_errors = classify_discards_and_errors(
            actual_obs=actual_obs,
            predictions=predictions,
            postselected_observables_mask=self.task.postselected_observables_mask,
            out_count_observable_error_combos=(
                custom_counts if self.count_observable_error_combos else None
            ),
            num_obs=self.num_obs,
        )

        custom_counts["low_confidence_shots"] += (
            self.compiled_decoder.last_batch_low_confidence_shots
        )

        return AnonTaskStats(
            shots=num_shots,
            errors=num_errors,
            discards=num_discards_1 + num_discards_2,
            seconds=time.monotonic() - t0,
            custom_counts=custom_counts,
        )


def get_tesseract_decoder(beam_setup: str) -> Dict[str, TesseractSampler]:
    if beam_setup == "long":
        kwargs = TESSERACT_LONG_KWARGS
    elif beam_setup == "short":
        kwargs = TESSERACT_SHORT_KWARGS
    else:
        raise ValueError(f"unknown beam_setup: {beam_setup!r}")

    decoder = TesseractSinterDecoder(**kwargs)
    sampler = TesseractSampler(
        decoder=decoder,
        count_detection_events=False,
        count_observable_error_combos=False,
        tmp_dir=None,
    )
    return {f"tesseract_{beam_setup}": sampler}
