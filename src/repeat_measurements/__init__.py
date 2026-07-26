"""Repeated-measurement Floquet circuit transform.

Repeats each mid-circuit measurement layer r times and inserts comparison
detectors between consecutive copies. Original DETECTOR and OBSERVABLE_INCLUDE
rec references are remapped so the resulting circuit preserves detector
semantics and (optionally) matchability.
"""

from repeat_measurements.repeat import repeat_measurements

__all__ = ["repeat_measurements"]
