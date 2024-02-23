"""The Cosmology API compatibility wrapper for CAMB."""

from __future__ import annotations

__all__ = ["StandardCosmologyWrapper"]

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy import vectorize

from . import constants
from ._core import CosmologyWrapper, InputT, NDFloating

_MPCS_KM_TO_GYR = np.array("978.5", dtype=np.float64)  # [Mpc s / km -> Gyr]


@dataclass(frozen=True)
class StandardCosmologyWrapper(CosmologyWrapper):
    """FLRW Cosmology API wrapper for CAMB cosmologies."""

    def __post_init__(self) -> None:
        """Run-time post-processing.

        Note that if this module is c-compiled (e.g. with :mod:`mypyc`) that
        the type of ``self.cosmo`` must be ``CAMBdata`` at object creation
        and cannot be later processed here.
        """
        super().__post_init__()

        self._cosmo_fn: dict[str, Any]
        object.__setattr__(
            self,
            "_cosmo_fn",
            {
                "comoving_distance": vectorize(
                    self.cosmo.comoving_radial_distance,
                    excluded={"tol"},
                ),
                "angular_diameter_distance": vectorize(
                    self.cosmo.angular_diameter_distance
                ),
                "luminosity_distance": vectorize(self.cosmo.luminosity_distance),
                "hubble_parameter": vectorize(self.cosmo.hubble_parameter),
                "get_Omega": vectorize(self.cosmo.get_Omega, excluded={"var"}),
            },
        )

    # ----------------------------------------------
    # TotalComponent

    @property
    def Omega_tot0(self) -> NDFloating:
        r"""Omega total; the total density/critical density at z=0.

        .. math::

            \Omega_{\rm tot} = \Omega_{\rm m} + \Omega_{\rm \gamma} +
            \Omega_{\rm \nu} + \Omega_{\rm de} + \Omega_{\rm k}
        """
        return np.array(self.cosmo.get_Omega("tot", z=0))

    def Omega_tot(self, z: InputT, /) -> NDFloating:
        r"""Redshift-dependent total density parameter.

        This is the sum of the matter, radiation, neutrino, dark energy, and
        curvature density parameters.

        .. math::

            \Omega_{\rm tot} = \Omega_{\rm m} + \Omega_{\rm \gamma} +
            \Omega_{\rm \nu} + \Omega_{\rm de} + \Omega_{\rm k}
        """
        return np.asarray(self._cosmo_fn["get_Omega"]("tot", z))

    # ----------------------------------------------
    # CurvatureComponent

    @property
    def Omega_k0(self) -> NDFloating:
        """Omega curvature; the effective curvature density/critical density at z=0."""
        return np.asarray(self._params.omk)

    def Omega_k(self, z: InputT, /) -> NDFloating:
        """Redshift-dependent curvature density parameter."""
        return np.asarray(self._cosmo_fn["get_Omega"]("K", z))

    # ----------------------------------------------
    # MatterComponent

    @property
    def Omega_m0(self) -> NDFloating:
        """Matter density at z=0."""
        return np.asarray(self._params.omegam)

    def Omega_m(self, z: InputT, /) -> NDFloating:
        """Redshift-dependent non-relativistic matter density parameter.

        Notes
        -----
        This does not include neutrinos, even if non-relativistic at the
        redshift of interest; see `Onu`.

        """
        return np.asarray(
            self._cosmo_fn["get_Omega"]("cdm", z)
            + self._cosmo_fn["get_Omega"]("baryon", z),
        )

    # ----------------------------------------------
    # BaryonComponent

    @property
    def Omega_b0(self) -> NDFloating:
        """Baryon density at z=0."""
        return np.asarray(self._params.omegab)

    def Omega_b(self, z: InputT, /) -> NDFloating:
        """Redshift-dependent baryon density parameter.

        Raises
        ------
        ValueError
            If ``Ob0`` is `None`.

        """
        return np.asarray(self._cosmo_fn["get_Omega"]("baryon", z))

    # ----------------------------------------------
    # NeutrinoComponent

    @property
    def Omega_nu0(self) -> NDFloating:
        """Omega nu; the density/critical density of neutrinos at z=0."""
        return np.asarray(
            self.cosmo.get_Omega("neutrino", z=0) + self.cosmo.get_Omega("nu", z=0)
        )

    @property
    def Neff(self) -> NDFloating:
        """Effective number of neutrino species."""
        return np.asarray(self._params.N_eff)

    @property
    def m_nu(self) -> tuple[NDFloating, ...]:
        """Neutrino mass in eV."""
        return tuple(np.array(self.cosmo.nu_masses))

    def Omega_nu(self, z: InputT, /) -> NDFloating:
        r"""Redshift-dependent neutrino density parameter."""
        return np.asarray(
            self._cosmo_fn["get_Omega"]("neutrino", z)
            + self._cosmo_fn["get_Omega"]("nu", z),
        )

    # ----------------------------------------------
    # DarkEnergyComponent

    @property
    def Omega_de0(self) -> NDFloating:
        """Dark energy density at z=0."""
        return np.asarray(self.cosmo.omega_de)

    def Omega_de(self, z: InputT, /) -> NDFloating:
        """Redshift-dependent dark energy density parameter."""
        return np.asarray(self._cosmo_fn["get_Omega"]("de", z))

    # ----------------------------------------------
    # DarkMatterComponent

    @property
    def Omega_dm0(self) -> NDFloating:
        """Omega dark matter; dark matter density/critical density at z=0."""
        return np.asarray(self._params.omegac)

    def Omega_dm(self, z: InputT, /) -> NDFloating:
        """Redshift-dependent dark matter density parameter.

        Notes
        -----
        This does not include neutrinos, even if non-relativistic at the
        redshift of interest.

        """
        return np.asarray(self._cosmo_fn["get_Omega"]("cdm", z))

    # ----------------------------------------------
    # PhotonComponent

    @property
    def Omega_gamma0(self) -> NDFloating:
        """Omega gamma; the density/critical density of photons at z=0."""
        return np.asarray(self.cosmo.get_Omega("photon", z=0))

    def Omega_gamma(self, z: InputT, /) -> NDFloating:
        """Redshift-dependent photon density parameter."""
        return np.asarray(self._cosmo_fn["get_Omega"]("photon", z))

    # ----------------------------------------------
    # CriticalDensity

    @property
    def critical_density0(self) -> NDFloating:
        """Critical density at z = 0 in Msol Mpc-3."""
        # H0 is in (km/s) Mpc-1; G is in pc Msol-1 (km/s)^2
        # so we pick up a factor of 1e6 to get Msol Mpc-3
        return np.array(3e6 * self.H0**2 / (8 * np.pi * constants.G))

    def critical_density(self, z: InputT, /) -> NDFloating:
        """Redshift-dependent critical density in Msol Mpc-3."""
        return np.array(3e6 * self.H(z) ** 2 / (8 * np.pi * constants.G))

    # ----------------------------------------------
    # HubbleParameter

    @property
    def H0(self) -> NDFloating:
        """Hubble constant at z=0 in km s-1 Mpc-1."""
        return np.array(self._params.H0)

    @property
    def hubble_distance(self) -> NDFloating:
        """Hubble distance in Mpc."""
        return np.array(constants.c / self.H0)

    @property
    def hubble_time(self) -> NDFloating:
        """Hubble time in Gyr."""
        return np.array(_MPCS_KM_TO_GYR / self.H0)

    def H(self, z: InputT, /) -> NDFloating:
        """Hubble function :math:`H(z)` in km s-1 Mpc-1."""  # noqa: D402
        return np.array(self._cosmo_fn["hubble_parameter"](z))

    def H_over_H0(self, z: InputT, /) -> NDFloating:
        """Standardised Hubble function :math:`E(z) = H(z)/H_0`."""
        return self.H(z) / self.H0

    # ----------------------------------------------
    # Scale factor

    @property
    def scale_factor0(self) -> NDFloating:
        """Scale factor at z=0."""
        return np.asarray(1.0)

    def scale_factor(self, z: InputT, /) -> NDFloating:
        """Redshift-dependenct scale factor :math:`a = a_0 / (1 + z)`."""
        return np.asarray(self.scale_factor0 / (z + 1))

    # ----------------------------------------------
    # Temperature

    @property
    def T_cmb0(self) -> NDFloating:
        """Temperature of the CMB at z=0."""
        return np.asarray(self._params.TCMB)

    def T_cmb(self, z: InputT, /) -> NDFloating:
        """Temperature of the CMB at redshift ``z``."""
        return self.T_cmb0 * (z + 1)

    # ----------------------------------------------
    # Time

    def age(self, z: InputT, /) -> NDFloating:
        """Age of the universe in Gyr at redshift ``z``."""
        return np.asarray(self.cosmo.physical_time(z))

    def lookback_time(self, z: InputT, /) -> NDFloating:
        """Lookback time to redshift ``z`` in Gyr.

        The lookback time is the difference between the age of the Universe now
        and the age at redshift ``z``.
        """
        # TODO: cache the age of the universe
        return np.asarray(self.cosmo.physical_time(0) - self.cosmo.physical_time(z))

    # ----------------------------------------------
    # Comoving distance

    def comoving_distance(self, z: InputT, zp: InputT | None = None, /) -> NDFloating:
        r"""Comoving line-of-sight distance :math:`d_c(z1, z2)` in Mpc.

        The comoving distance along the line-of-sight between two objects
        remains constant with time for objects in the Hubble flow.

        Parameters
        ----------
        z, zp : ndarray, positional-only
            Input redshifts. If ``zp`` is `None` (default), then the distance
            :math:`d_c(0, z)` is returned, otherwise the distance :math:`d_c(z,
            zp)` is returned.

        Returns
        -------
        ndarray
            The comoving distance :math:`d_c(z1, z2)` in Mpc, where ``(z1, z2)``
            is (0, `z`) if `zp` is `None` else (`z`, `zp`).

        """
        # TODO: have a way to set the tolerance
        z1, z2 = (0, z) if zp is None else (z, zp)
        return np.asarray(
            self._cosmo_fn["comoving_distance"](z2, tol=0.0001)
        ) - np.asarray(self._cosmo_fn["comoving_distance"](z1, tol=0.0001))

    def comoving_transverse_distance(
        self, z: InputT, zp: InputT | None = None, /
    ) -> NDFloating:
        r"""Transverse comoving distance :math:`d_M(z1, z2)` in Mpc.

        This value is the transverse comoving distance at redshift ``z``
        corresponding to an angular separation of 1 radian. This is the same as
        the comoving distance if :math:`\Omega_k` is zero (as in the current
        concordance Lambda-CDM model).

        Parameters
        ----------
        z, zp : ndarray, positional-only
            Input redshifts. If ``zp`` is `None` (default), then the distance
            :math:`d_M(0, z)` is returned, otherwise the distance :math:`d_M(z,
            zp)` is returned.

        Returns
        -------
        ndarray
            The comoving transverse distance :math:`d_M(z1, z2)` in Mpc, where
            ``(z1, z2)`` is (0, `z`) if `zp` is `None` else (`z`, `zp`).

        """
        z1, z2 = (0, z) if zp is None else (z, zp)
        return (
            self._cosmo_fn["angular_diameter_distance"](z2)
            - self._cosmo_fn["angular_diameter_distance"](z1)
        ) * (z2 + 1)

    def _comoving_volume_flat(self, z: InputT, /) -> NDFloating:
        """Comoving volume in cubic Mpc for flat cosmologies."""
        return 4.0 / 3.0 * np.pi * self.comoving_distance(z) ** 3

    def _comoving_volume_positive(self, z: InputT, /) -> NDFloating:
        dh = self.hubble_distance
        x = self.comoving_transverse_distance(z) / dh
        term1 = 4.0 * np.pi * dh**3 / (2.0 * self.Omega_k0)
        term2 = x * np.sqrt(1 + self.Omega_k0 * (x) ** 2)
        term3 = np.sqrt(np.abs(self.Omega_k0)) * x

        return term1 * (
            term2 - 1.0 / np.sqrt(np.abs(self.Omega_k0)) * np.arcsinh(term3)
        )

    def _comoving_volume_negative(self, z: InputT, /) -> NDFloating:
        dh = self.hubble_distance
        x = self.comoving_transverse_distance(z) / dh
        term1 = 4.0 * np.pi * dh**3 / (2.0 * self.Omega_k0)
        term2 = x * np.sqrt(1 + self.Omega_k0 * (x) ** 2)
        term3 = np.sqrt(np.abs(self.Omega_k0)) * x

        return term1 * (term2 - 1.0 / np.sqrt(np.abs(self.Omega_k0)) * np.arcsin(term3))

    def comoving_volume(self, z: InputT, zp: InputT | None = None, /) -> NDFloating:
        r"""Comoving volume in cubic Mpc.

        This is the volume of the universe encompassed by redshifts less than
        ``z``. For the case of :math:`\Omega_k = 0` it is a sphere of radius
        `comoving_distance` but it is less intuitive if :math:`\Omega_k` is not.

        Parameters
        ----------
        z, zp : ndarray, positional-only
            Input redshifts. If ``zp`` is `None` (default), then the
            volume :math:`V_c(0, z)` is returned, otherwise the
            volume :math:`V_c(z, zp)` is returned.

        Returns
        -------
        ndarray
            The comoving volume :math:`V_c(z1, z2)` in Mpc, where
            ``(z1, z2)`` is (0, `z`) if `zp` is `None` else (`z`, `zp`).

        """
        z1, z2 = (0, z) if zp is None else (z, zp)
        if self.Omega_k0 == 0:
            vc = self._comoving_volume_flat(z2) - self._comoving_volume_flat(z1)
        elif self.Omega_k0 > 0:
            vc = self._comoving_volume_positive(z2) - self._comoving_volume_positive(z1)
        else:
            vc = self._comoving_volume_negative(z2) - self._comoving_volume_negative(z1)
        return vc

    def differential_comoving_volume(
        self, z: InputT, zp: InputT | None = None, /
    ) -> NDFloating:
        r"""Differential comoving volume in cubic Mpc per steradian.

        If :math:`V_c` is the comoving volume of a redshift slice with solid
        angle :math:`\Omega`, this function returns

        .. math::

            \mathtt{dvc(z)}
            = \frac{1}{d_H^3} \, \frac{dV_c}{d\Omega \, dz}
            = \frac{x_M^2(z)}{E(z)}
            = \frac{\mathtt{xm(z)^2}}{\mathtt{ef(z)}} \;.

        Parameters
        ----------
        z, zp : ndarray, positional-only
            Input redshifts. If ``zp`` is `None` (default), then the
            differential volume :math:`dV_c(0, z)` is returned, otherwise the
            differential volume :math:`dV_c(z, zp)` is returned.

        Returns
        -------
        ndarray
            The differential comoving volume :math:`dV_c(z1, z2)` in Mpc,
            where ``(z1, z2)`` is (0, `z`) if `zp` is `None` else (`z`, `zp`).

        """
        z1, z2 = (0, zp) if zp is None else (z, zp)
        return (
            self.comoving_transverse_distance(z2) / self.hubble_distance
        ) ** 2 / self.H_over_H0(z2) - (
            self.comoving_transverse_distance(z1) / self.hubble_distance
        ) ** 2 / self.H_over_H0(z1)

    # ----------------------------------------------
    # Angular diameter distance

    def angular_diameter_distance(
        self, z: InputT, zp: InputT | None = None, /
    ) -> NDFloating:
        """Angular diameter distance :math:`d_A(z)` in Mpc.

        This gives the proper (sometimes called 'physical') transverse distance
        corresponding to an angle of 1 radian for an object at redshift ``z``
        ([1]_, [2]_, [3]_).

        Parameters
        ----------
        z, zp : ndarray, positional-only
            Input redshifts. If ``zp`` is `None` (default), then the distance
            :math:`d_A(0, z)` is returned, otherwise the distance :math:`d_A(z,
            zp)` is returned.

        Returns
        -------
        ndarray
            The angular diameter distance :math:`d_A(z1, z2)` in Mpc, where
            ``(z1, z2)`` is (0, `z`) if `zp` is `None` else (`z`, `zp`).

        References
        ----------
        .. [1] Weinberg, 1972, pp 420-424; Weedman, 1986, pp 421-424.
        .. [2] Weedman, D. (1986). Quasar astronomy, pp 65-67.
        .. [3] Peebles, P. (1993). Principles of Physical Cosmology, pp 325-327.

        """
        z1, z2 = (0, z) if zp is None else (z, zp)
        return np.asarray(
            self._cosmo_fn["angular_diameter_distance"](z2)
            - self._cosmo_fn["angular_diameter_distance"](z1)
        )

    # ----------------------------------------------
    # Luminosity distance

    def luminosity_distance(self, z: InputT, zp: InputT | None = None, /) -> NDFloating:
        """Redshift-dependent luminosity distance :math:`d_L(z1, z2)` in Mpc.

        This is the distance to use when converting between the bolometric flux
        from an object at redshift ``z`` and its bolometric luminosity [1]_.

        Parameters
        ----------
        z, zp : ndarray, positional-only
            Input redshifts. If ``zp`` is `None` (default), then the
            distance :math:`d_L(0, z)` is returned, otherwise the
            distance :math:`d_L(z, zp)` is returned.

        Returns
        -------
        ndarray
            The luminosity distance :math:`d_L(z1, z2)` in Mpc, where
            ``(z1, z2)`` is (0, `z`) if `zp` is `None` else (`z`, `zp`).

        References
        ----------
        .. [1] Weinberg, 1972, pp 420-424; Weedman, 1986, pp 60-62.

        """
        z1, z2 = (0, z) if zp is None else (z, zp)
        return self._cosmo_fn["luminosity_distance"](z2) - self._cosmo_fn[
            "luminosity_distance"
        ](z1)
