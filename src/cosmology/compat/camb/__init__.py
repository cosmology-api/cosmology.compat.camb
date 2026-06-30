"""Cosmology API compatibility layer for CAMB."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

TYPE_CHECKING = False
if TYPE_CHECKING:
    from types import ModuleType
    from typing import TypeAlias

    import jaxtyping
    from array_api_strict._array_object import Array
    from numpy.typing import NDArray

    from camb import CAMBdata, CAMBparams

    FloatArray: TypeAlias = NDArray[np.float64] | jaxtyping.Array | Array


@dataclass(frozen=True, slots=True)
class Cosmology:
    """Cosmology API wrapper for CAMB *pars* and *results*."""

    data: CAMBdata
    params: CAMBparams = field(init=False)
    xp: ModuleType = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "params", self.data.Params)
        object.__setattr__(self, "xp", np)

    def set_xp(self, xp: ModuleType) -> None:
        """Setter for the array backend (`xp`), defaults to NumPy."""
        object.__setattr__(self, "xp", xp)

    @property
    def h(self) -> FloatArray:
        """Little h."""
        return self.xp.asarray(self.params.h)

    @property
    def H0(self) -> FloatArray:
        """Hubble constant."""
        return self.xp.asarray(self.params.H0)

    @property
    def Omega_m0(self) -> FloatArray:
        """Total matter today, excluding massive neutrinos."""
        return self.xp.asarray(self.params.omegam)

    @property
    def Omega_de0(self) -> FloatArray:
        """Dark energy today."""
        return self.xp.asarray(self.data.omega_de)

    @property
    def Omega_k0(self) -> FloatArray:
        """Curvature today."""
        return self.xp.asarray(self.params.omk)

    @property
    def hubble_distance(self) -> FloatArray:
        """Hubble distance."""
        return self.xp.asarray(299792.458 / self.params.H0)

    @property
    def critical_density0(self) -> FloatArray:
        """Critical density today in Msol Mpc-3."""
        # gravitational constant kappa = 8pi G/c^2 in Mpc Msol-1
        # uses nominal value of (G Msol) following IAU 2015
        kappa = 1.202706180375887e-18
        return self.xp.asarray(self.data.grhocrit / kappa)

    def H(self, z: FloatArray | float) -> FloatArray:
        """Hubble parameter at redshift *z*."""
        return self.xp.asarray(self.data.hubble_parameter(np.asarray(z)))

    def Omega_m(self, z: FloatArray | float) -> FloatArray:
        """Total matter, excluding massive neutrinos, at redshift *z*."""
        z_np = np.asarray(z)
        return self.xp.asarray(
            self.data.get_Omega("baryon", z_np)
            + self.data.get_Omega("cdm", z_np)
            + self.data.get_Omega("nu", z_np)
        )

    def Omega_de(self, z: FloatArray | float) -> FloatArray:
        """Dark energy at redshift *z*."""
        return self.xp.asarray(self.data.get_Omega("de", np.asarray(z)))

    def Omega_k(self, z: FloatArray | float) -> FloatArray:
        """Curvature at redshift *z*."""
        return self.xp.asarray(self.data.get_Omega("K", np.asarray(z)))

    def comoving_distance(
        self,
        z: FloatArray | float,
        z2: FloatArray | float | None = None,
    ) -> FloatArray:
        """Comoving distance at redshift *z*.

        If *z2* is given, computes the comoving distance between
        redshifts *z* and *z2*.

        """
        z_np = np.asarray(z)
        if z2 is not None:
            return self.xp.asarray(
                self.data.comoving_radial_distance(np.asarray(z2))
                - self.data.comoving_radial_distance(z_np),
            )
        return self.xp.asarray(self.data.comoving_radial_distance(z_np))

    def inv_comoving_distance(self, x: FloatArray | float) -> FloatArray:
        """Return redshift at which the comoving distance is *x*."""
        return self.xp.asarray(
            self.data.redshift_at_comoving_radial_distance(np.asarray(x))
        )

    def angular_diameter_distance(
        self,
        z: FloatArray | float,
        z2: FloatArray | float | None = None,
    ) -> FloatArray:
        """Angular diameter distance at redshift *z*.

        If *z2* is given, computes the Angular diameter distance between
        redshifts *z* and *z2*.

        """
        z_np = np.asarray(z)
        if z2 is not None:
            return self.xp.asarray(
                self.data.angular_diameter_distance2(z_np, np.asarray(z2))
            )
        return self.xp.asarray(self.data.angular_diameter_distance(z_np))

    def H_over_H0(self, z: FloatArray | float) -> FloatArray:
        """Standardised Hubble function :math:`E(z) = H(z)/H_0`."""
        return self.H(np.asarray(z)) / self.xp.asarray(self.H0)

    def transverse_comoving_distance(
        self,
        z: FloatArray | float,
        z2: FloatArray | float | None = None,
    ) -> FloatArray:
        """Transverse comoving distance at redshift *z*.

        If *z2* is given, computes the transverse comoving distance between
        redshifts *z* and *z2*.

        """
        z_np = np.asarray(z)
        if z2 is not None:
            z2_np = np.asarray(z2)
            return self.xp.asarray(
                (1 + z2_np)
                * self.data.angular_diameter_distance2(z_np, np.asarray(z2_np))
            )
        return self.xp.asarray((1 + z_np) * self.data.angular_diameter_distance(z_np))
