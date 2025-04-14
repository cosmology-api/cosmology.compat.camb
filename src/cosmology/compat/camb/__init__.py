"""Cosmology API compatibility layer for CAMB."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

TYPE_CHECKING = False
if TYPE_CHECKING:
    from typing import TypeAlias

    from numpy.typing import NDArray

    from camb import CAMBdata, CAMBparams

    Array: TypeAlias = NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class Cosmology:
    """Cosmology API wrapper for CAMB *pars* and *results*."""

    data: CAMBdata
    params: CAMBparams = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "params", self.data.Params)

    @property
    def h(self) -> Array:
        """Little h."""
        return np.array(self.params.h)

    @property
    def H0(self) -> Array:
        """Hubble constant."""
        return np.array(self.params.H0)

    @property
    def Omega_m0(self) -> Array:
        """Total matter today, excluding massive neutrinos."""
        return np.array(self.params.omegam)

    @property
    def Omega_de0(self) -> Array:
        """Dark energy today."""
        return np.array(self.data.omega_de)

    @property
    def Omega_k0(self) -> Array:
        """Curvature today."""
        return np.array(self.params.omk)

    @property
    def hubble_distance(self) -> Array:
        """Hubble distance."""
        return np.array(299792.458 / self.params.H0)

    @property
    def critical_density0(self) -> Array:
        """Critical density today in Msol Mpc-3."""
        # gravitational constant kappa = 8pi G/c^2 in Mpc Msol-1
        # uses nominal value of (G Msol) following IAU 2015
        kappa = 1.202706180375887e-18
        return np.array(self.data.grhocrit / kappa)

    def H(self, z: Array | float) -> Array:
        """Hubble parameter at redshift *z*."""
        return np.array(self.data.hubble_parameter(z))

    def Omega_m(self, z: Array | float) -> Array:
        """Total matter, excluding massive neutrinos, at redshift *z*."""
        return np.array(
            self.data.get_Omega("baryon", z)
            + self.data.get_Omega("cdm", z)
            + self.data.get_Omega("nu", z)
        )

    def Omega_de(self, z: Array | float) -> Array:
        """Dark energy at redshift *z*."""
        return np.array(self.data.get_Omega("de", z))

    def Omega_k(self, z: Array | float) -> Array:
        """Curvature at redshift *z*."""
        return np.array(self.data.get_Omega("K", z))

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
                self.data.comoving_radial_distance(z2)
                - self.data.comoving_radial_distance(z),
            )
        return np.array(self.data.comoving_radial_distance(z))

    def inv_comoving_distance(self, x: Array | float) -> Array:
        """Return redshift at which the comoving distance is *x*."""
        return np.array(self.data.redshift_at_comoving_radial_distance(x))

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
            return np.array(self.data.angular_diameter_distance2(z, z2))
        return np.array(self.data.angular_diameter_distance(z))

    def H_over_H0(self, z: Array | float) -> Array:
        """Standardised Hubble function :math:`E(z) = H(z)/H_0`."""
        return self.H(z) / self.H0

    def transverse_comoving_distance(
        self,
        z: Array | float,
        z2: Array | float | None = None,
    ) -> Array:
        """Transverse comoving distance at redshift *z*.

        If *z2* is given, computes the transverse comoving distance between
        redshifts *z* and *z2*.

        """
        if z2 is not None:
            return np.array((1 + z2) * self.data.angular_diameter_distance2(z, z2))
        return np.array((1 + z) * self.data.angular_diameter_distance(z))
