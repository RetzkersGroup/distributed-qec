# Threshold Comparison: flexible-pairwise noise model

## Goal

Compare thresholds for different noise regimes on native pairwise-measurement circuits (no Bell expansion). The noise model separates three independent error sources to study which channel dominates the threshold.

## Noise Model: `flexible_pairwise`

Applied directly to circuits containing `MPP` (2-body), `R`/`RX`, and `M`/`MX` instructions using `_cultiv_noise.NoiseModel.noisy_circuit()`.

### Parameters

| Parameter | Channel | Applied to |
|-----------|---------|-----------|
| `p_m_bulk` | Measurement flip (`flip_result`) | 2-body pairwise measurements (XX, YY, ZZ) |
| `p_spam` | Reset error + single-qubit measurement flip | `X_ERROR` after R, `Z_ERROR` after RX, `flip_result` on M/MX |
| `p_depolarizing` | `DEPOLARIZE1` after measurement | 2-body pairwise measurements (XX, YY, ZZ) |

No idle noise, no single-qubit gate noise (H/S are treated as noiseless logical basis changes).

We use `DEPOLARIZE1` because we assume the underlying implementation is Bell-measurement which do not introduce correlated errors.

### Usage

```python
from _cultiv_noise import NoiseModel

model = NoiseModel.flexible_pairwise(p_m_bulk=0.01, p_spam=0.005, p_depolarizing=0.02)
noisy_circuit = model.noisy_circuit(raw_circuit.without_noise())
```

Or via the dispatcher with named profiles:

```python
from noise_utils import get_noise_model

model = get_noise_model("flexible_pairwise", p=0.01, profile="measurement-limited")
```

## Comparison Method: Named Ratio Profiles

Each profile fixes the ratio among the 3 parameters, then we sweep the overall scale `p`:

| Profile | p_m_bulk | p_spam | p_depolarizing |
|---------|----------|--------|----------------|
| `balanced` | p | p | p |
| `measurement-limited` | p | 0 | 0 |
| `depolarizing-limited` | 0 | 0 | p |
| `spam-limited` | 0 | p | 0 |

The threshold for each profile tells us which noise channel is most damaging to the code.

## Pipeline

```
input circuit (.stim, noiseless)
    → NoiseModel.flexible_pairwise(...).noisy_circuit(raw)
    → sinter.Task (in-memory)
    → sinter.collect() with pymatching
    → plot: logical error rate vs p, grouped by (distance, profile)
```

See `notebooks/threshold_comparison.ipynb` for the full workflow.
