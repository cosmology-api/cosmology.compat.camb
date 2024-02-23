"""Test the Cosmology API compat library."""

from __future__ import annotations

from typing import Any

from hypothesis.extra import numpy as npst


# Hypothesis strategy for generating arrays of redshifts
def z_arr_st(
    *,
    allow_nan: bool = False,
    min_value: float | None = 0,
    max_value: float | None = None,
    **kwargs: Any,
) -> Any:
    """Hypothesis strategy for generating arrays of redshifts."""
    return npst.arrays(
        # TODO: do we want to test float16?
        npst.floating_dtypes(sizes=(32, 64)),
        npst.array_shapes(min_dims=1),
        elements={
            "allow_nan": allow_nan,
            "min_value": min_value,
            "max_value": max_value,
        }
        | kwargs,
    )
