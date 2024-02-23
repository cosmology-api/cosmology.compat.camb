"""Test the Cosmology API compat library."""

from __future__ import annotations

import numpy as np
from hypothesis import given

from cosmology.api import (
    BaryonComponent,
    CurvatureComponent,
    DarkEnergyComponent,
    DarkMatterComponent,
    MatterComponent,
    NeutrinoComponent,
    PhotonComponent,
    TotalComponent,
)

from .conftest import z_arr_st

################################################################################
# TESTS
################################################################################


class TotalComponent_Test:
    def test_wrapper_is_compliant(self, wrapper):
        """Test that AstropyCosmology is a BackgroundCosmologyWrapper."""
        if hasattr(super(), "test_wrapper_is_compliant"):
            super().test_wrapper_is_compliant(wrapper)

        assert isinstance(wrapper, TotalComponent)

    def test_Omega_tot0(self, wrapper, cosmo):
        """Test that the wrapper has the same Otot0 as the wrapped object."""
        assert wrapper.Omega_tot0 == cosmo.get_Omega("tot", z=0)
        assert isinstance(wrapper.Omega_tot0, np.ndarray)

    @given(z_arr_st())
    def test_Otot(self, wrapper, vcosmo, z):
        """Test that the wrapper's Otot is the same as the wrapped object's."""
        omega = wrapper.Omega_tot(z)
        assert np.array_equal(omega, vcosmo.get_Omega("tot", z))
        assert isinstance(omega, np.ndarray)


class CurvatureComponent_Test:
    def test_wrapper_is_compliant(self, wrapper):
        """Test that AstropyCosmology is a BackgroundCosmologyWrapper."""
        if hasattr(super(), "test_wrapper_is_compliant"):
            super().test_wrapper_is_compliant(wrapper)

        assert isinstance(wrapper, CurvatureComponent)

    def test_Omega_k0(self, wrapper, cosmo, params):
        """Test that the wrapper has the same Omega_k0 as the wrapped object."""
        assert wrapper.Omega_k0 == params.omk
        assert np.allclose(wrapper.Omega_k0, cosmo.get_Omega("K", z=0))
        assert isinstance(wrapper.Omega_k0, np.ndarray)

    @given(z_arr_st())
    def test_Omega_k(self, wrapper, vcosmo, z):
        """Test that the wrapper's Omega_k is the same as the wrapped object's."""
        omega = wrapper.Omega_k(z)
        assert np.array_equal(omega, vcosmo.get_Omega("K", z))
        assert isinstance(omega, np.ndarray)


class MatterComponent_Test:
    def test_wrapper_is_compliant(self, wrapper):
        """Test that AstropyCosmology is a BackgroundCosmologyWrapper."""
        if hasattr(super(), "test_wrapper_is_compliant"):
            super().test_wrapper_is_compliant(wrapper)

        assert isinstance(wrapper, MatterComponent)

    def test_Omega_m0(self, wrapper, cosmo, params):
        """Test that the wrapper has the same Om0 as the wrapped object."""
        assert wrapper.Omega_m0 == params.omegam
        # assert np.allclose(
        #     wrapper.Omega_m0,
        # )  # TODO: why is this different than the above?
        assert isinstance(wrapper.Omega_m0, np.ndarray)

    # TODO: why do these fail for z-> inf?
    @given(z_arr_st(max_value=1e9))
    def test_Omega_m(self, wrapper, vcosmo, z):
        """Test that the wrapper's Om is the same as the wrapped object's."""
        omega = wrapper.Omega_m(z)
        assert np.array_equal(
            omega, vcosmo.get_Omega("cdm", z) + vcosmo.get_Omega("baryon", z)
        )
        assert isinstance(omega, np.ndarray)


class BaryonComponent_Test:
    def test_wrapper_is_compliant(self, wrapper):
        """Test that AstropyCosmology is a BackgroundCosmologyWrapper."""
        if hasattr(super(), "test_wrapper_is_compliant"):
            super().test_wrapper_is_compliant(wrapper)

        assert isinstance(wrapper, BaryonComponent)

    def test_Omega_b0(self, wrapper, cosmo, params):
        """Test that the wrapper has the same Omega_b0 as the wrapped object."""
        assert wrapper.Omega_b0 == params.omegab
        assert np.allclose(wrapper.Omega_b0, cosmo.get_Omega("baryon", z=0))
        assert isinstance(wrapper.Omega_b0, np.ndarray)

    @given(z_arr_st(max_value=1e9))
    def test_Omega_b(self, wrapper, vcosmo, z):
        """Test that the wrapper's Omega_b is the same as the wrapped object's."""
        omega = wrapper.Omega_b(z)
        assert np.array_equal(omega, vcosmo.get_Omega("baryon", z))
        assert isinstance(omega, np.ndarray)


class NeutrinoComponent_Test:
    def test_wrapper_is_compliant(self, wrapper):
        """Test that AstropyCosmology is a BackgroundCosmologyWrapper."""
        if hasattr(super(), "test_wrapper_is_compliant"):
            super().test_wrapper_is_compliant(wrapper)

        assert isinstance(wrapper, NeutrinoComponent)

    def test_Omega_nu0(self, wrapper, cosmo):
        """Test that the wrapper has the same Omega_nu0 as the wrapped object."""
        assert wrapper.Omega_nu0 == cosmo.get_Omega("neutrino", z=0) + cosmo.get_Omega(
            "nu", z=0
        )
        assert isinstance(wrapper.Omega_nu0, np.ndarray)

    def test_Neff(self, wrapper, params):
        """Test that the wrapper has the same Neff as the wrapped object."""
        assert wrapper.Neff == params.N_eff
        assert isinstance(wrapper.Neff, np.ndarray)

    def test_m_nu(self, wrapper, cosmo):
        """Test that the wrapper has the same m_nu as the wrapped object."""
        assert all(np.equal(w, c) for w, c in zip(wrapper.m_nu, cosmo.nu_masses))
        # assert isinstance(wrapper.m_nu, tuple) for m in wrapper.m_nu)

    @given(z_arr_st(max_value=1e9))
    def test_Omega_nu(self, wrapper, vcosmo, z):
        """Test that the wrapper's Omega_nu is the same as the wrapped object's."""
        omega = wrapper.Omega_nu(z)
        assert np.array_equal(
            omega, vcosmo.get_Omega("neutrino", z) + vcosmo.get_Omega("nu", z)
        )
        assert isinstance(omega, np.ndarray)


class DarkEnergyComponent_Test:
    def test_wrapper_is_compliant(self, wrapper):
        """Test that AstropyCosmology is a BackgroundCosmologyWrapper."""
        if hasattr(super(), "test_wrapper_is_compliant"):
            super().test_wrapper_is_compliant(wrapper)

        assert isinstance(wrapper, DarkEnergyComponent)

    def test_Omega_de0(self, wrapper, cosmo):
        """Test that the wrapper has the same Omega_de0 as the wrapped object."""
        assert wrapper.Omega_de0 == cosmo.omega_de
        assert isinstance(wrapper.Omega_de0, np.ndarray)

    @given(z_arr_st())
    def test_Omega_de(self, wrapper, vcosmo, z):
        """Test that the wrapper's Omega_de is the same as the wrapped object's."""
        omega = wrapper.Omega_de(z)
        assert np.array_equal(omega, vcosmo.get_Omega("de", z))
        assert isinstance(omega, np.ndarray)


class DarkMatterComponent_Test:
    def test_wrapper_is_compliant(self, wrapper):
        """Test that AstropyCosmology is a BackgroundCosmologyWrapper."""
        if hasattr(super(), "test_wrapper_is_compliant"):
            super().test_wrapper_is_compliant(wrapper)

        assert isinstance(wrapper, DarkMatterComponent)

    def test_Omega_dm0(self, wrapper, cosmo, params):
        """Test that the wrapper has the same Omega_dm0 as the wrapped object."""
        assert wrapper.Omega_dm0 == params.omegac
        assert np.allclose(wrapper.Omega_dm0, cosmo.get_Omega("cdm", 0))
        assert isinstance(wrapper.Omega_dm0, np.ndarray)

    @given(z_arr_st(max_value=1e9))
    def test_Omega_dm(self, wrapper, vcosmo, z):
        """Test that the wrapper's Omega_dm is the same as the wrapped object's."""
        omega = wrapper.Omega_dm(z)
        assert np.array_equal(omega, vcosmo.get_Omega("cdm", z))
        assert isinstance(omega, np.ndarray)


class PhotonComponent_Test:
    def test_wrapper_is_compliant(self, wrapper):
        """Test that AstropyCosmology is a BackgroundCosmologyWrapper."""
        if hasattr(super(), "test_wrapper_is_compliant"):
            super().test_wrapper_is_compliant(wrapper)

        assert isinstance(wrapper, PhotonComponent)

    def test_Omega_gamma0(self, wrapper, cosmo):
        """Test that the wrapper has the same Omega_gamma0 as the wrapped object."""
        assert wrapper.Omega_gamma0 == cosmo.get_Omega("photon")
        assert isinstance(wrapper.Omega_gamma0, np.ndarray)

    @given(z_arr_st(max_value=1e9))
    def test_Omega_gamma(self, wrapper, vcosmo, z):
        """Test that the wrapper's Omega_gamma is the same as the wrapped object's."""
        omega = wrapper.Omega_gamma(z)
        assert np.array_equal(omega, vcosmo.get_Omega("photon", z))
        assert isinstance(omega, np.ndarray)
