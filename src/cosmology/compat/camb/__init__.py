"""Cosmology API compatibility layer for CAMB."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

TYPE_CHECKING = False
if TYPE_CHECKING:
    from typing import TypeAlias

    from numpy.typing import NDArray

    from camb import CAMBdata, CAMBparams

    Array: TypeAlias = NDArray[np.float64]


@dataclass
class Cosmology:
    """Cosmology API wrapper for CAMB *pars* and *results*."""

    pars: CAMBparams
    results: CAMBdata

    @property
    def h(self) -> Array:
        """Little h."""
        return np.array(self.pars.h)

    @property
    def H0(self) -> Array:
        """Hubble constant."""
        return np.array(self.pars.H0)

    @property
    def Omega_m0(self) -> Array:
        """Total matter today, excluding massive neutrinos."""
        return np.array(self.pars.omegam)

    @property
    def Omega_de0(self) -> Array:
        """Dark energy today."""
        return np.array(self.results.omega_de)

    @property
    def Omega_k0(self) -> Array:
        """Curvature today."""
        return np.array(self.pars.omk)

    @property
    def hubble_distance(self) -> Array:
        """Hubble distance."""
        return np.array(299792.458 / self.pars.H0)

    def H(self, z: Array | float) -> Array:
        """Hubble parameter at redshift *z*."""
        return np.array(self.results.hubble_parameter(z))

    def Omega_m(self, z: Array | float) -> Array:
        """Total matter, excluding massive neutrinos, at redshift *z*."""
        return np.array(
            self.results.get_Omega("baryon", z)
            + self.results.get_Omega("cdm", z)
            + self.results.get_Omega("nu", z)
        )

    def Omega_de(self, z: Array | float) -> Array:
        """Dark energy at redshift *z*."""
        return np.array(self.results.get_Omega("de", z))

    def Omega_k(self, z: Array | float) -> Array:
        """Curvature at redshift *z*."""
        return np.array(self.results.get_Omega("K", z))

    def comoving_distance(
        self,
        z: Array | float,
        z2: Array | float | None = None,
    ) -> Array:
        """Comoving distance at redshift *z*.

        If *z2* is given, computes the comoving distance between
        redshifts *z* and *z2*.

        """
        if z2 is not None:
            return np.array(
                self.results.comoving_radial_distance(z2)
                - self.results.comoving_radial_distance(z),
            )
        return np.array(self.results.comoving_radial_distance(z))

    def inv_comoving_distance(self, x: Array | float) -> Array:
        """Return redshift at which the comoving distance is *x*."""
        return np.array(self.results.redshift_at_comoving_radial_distance(x))

    def angular_diameter_distance(
        self,
        z: Array | float,
        z2: Array | float | None = None,
    ) -> Array:
        """Angular diameter distance at redshift *z*.

        If *z2* is given, computes the Angular diameter distance between
        redshifts *z* and *z2*.

        """
        if z2 is not None:
            return np.array(self.results.angular_diameter_distance2(z, z2))
        return np.array(self.results.angular_diameter_distance(z))

    def H_over_H0(self, z: Array | float) -> Array:
        """Standardised Hubble function :math:`E(z) = H(z)/H_0`."""
        return self.H(z) / self.H0
