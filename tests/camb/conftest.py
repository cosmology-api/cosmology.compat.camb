"""Fixtures shared across multiple test files."""

from __future__ import annotations

from typing import TYPE_CHECKING

import array_api_strict
import jax
import numpy as np
import pytest

if TYPE_CHECKING:
    from types import ModuleType

xp_available_backends: dict[str, ModuleType] = {
    "numpy": np,
    "array_api_strict": array_api_strict,
    "jax.numpy": jax.numpy,
}

# enable 64 bit numbers
jax.config.update("jax_enable_x64", val=True)
array_api_strict.set_array_api_strict_flags(api_version="2025.12")


@pytest.fixture(params=xp_available_backends.values(), scope="session")
def xp(request: pytest.FixtureRequest) -> ModuleType:
    """Fixture for array backend.

    Access array library functions using `xp.` in tests.

    """
    return request.param
