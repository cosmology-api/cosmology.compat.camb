"""The Cosmology API compatibility library for :mod:`camb`.

This library provides wrappers for CAMB cosmology objects to be compatible with
the Cosmology API. The available wrappers are:

- :class:`.StandardCosmology`: the Cosmology API wrapper for
  :mod:`camb`.


There are the following required objects for a Cosmology-API compatible library:

- constants: a module of constants. See :mod:`cosmology.compat.camb.constants`
  for details.
"""

from cosmology.compat.camb import constants
from cosmology.compat.camb._standard import StandardCosmologyWrapper

__all__ = [
    # Cosmology API
    "constants",
    # Wrappers
    "StandardCosmologyWrapper",
]
