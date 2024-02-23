"""CAMB cosmology constants.

From the :mod:`cosmology.api`, the list of required constants is:

- c: Speed of light in km s-1.
- G: Gravitational constant G in pc km2 s-2 Msol-1.
"""

import numpy as np
from camb.constants import c as c_ms

from cosmology.compat.camb._core import NDFloating

__all__ = ["c", "G"]

c: NDFloating = np.array(c_ms / 1e3)  # [km s-1]
# G: CODATA 2018 value
G: NDFloating = np.array(4.30091727003628e-3)  # [pc km2 s-2 Msol-1]
