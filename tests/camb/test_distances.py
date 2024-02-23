"""Test the Cosmology API compat library."""

from __future__ import annotations

import numpy as np
import pytest
from hypothesis import given, settings

from cosmology.api import DistanceMeasures

from .conftest import z_arr_st

################################################################################
# TESTS
################################################################################


class DistanceMeasures_Test:
    def test_wrapper_is_compliant(self, wrapper):
        """Test that AstropyCosmology is a BackgroundCosmologyWrapper."""
        if hasattr(super(), "test_wrapper_is_compliant"):
            super().test_wrapper_is_compliant(wrapper)

        assert isinstance(wrapper, DistanceMeasures)

    # =========================================================================

    def test_scale_factor0(self, wrapper):
        """Test that scale_factor(0) returns 1."""
        assert wrapper.scale_factor0 == 1
        assert isinstance(wrapper.scale_factor0, np.ndarray)

    @given(z_arr_st(min_value=None))
    def test_scale_factor(self, wrapper, cosmo, z):
        """Test that the wrapper's scale_factor is the same as the wrapped object's."""
        a = wrapper.scale_factor(z)
        assert np.array_equal(a, 1 / (1 + z))
        assert isinstance(a, np.ndarray)

    def test_T_cmb0(self, wrapper, params):
        """Test that the wrapper has the same Tcmb0 as the wrapped object."""
        assert wrapper.T_cmb0 == params.TCMB
        assert isinstance(wrapper.T_cmb0, np.ndarray)

    @given(z_arr_st())
    def test_T_cmb(self, wrapper, cosmo, params, z):
        """Test that the wrapper's Tcmb is the same as the wrapped object's."""
        T = wrapper.T_cmb(z)
        assert np.array_equal(T, params.TCMB * (1 + z))
        assert isinstance(T, np.ndarray)

    @given(z_arr_st())
    def test_age(self, wrapper, cosmo, z):
        """Test the wrapper's age."""
        age = wrapper.age(z)
        assert np.array_equal(age, cosmo.physical_time(z))
        assert isinstance(age, np.ndarray)

    @given(z_arr_st())
    def test_lookback_time(self, wrapper, cosmo, z):
        """Test the wrapper's lookback_time."""
        t = wrapper.lookback_time(z)
        assert np.array_equal(t, cosmo.physical_time(0) - cosmo.physical_time(z))
        assert isinstance(t, np.ndarray)

    @settings(deadline=500)
    @given(z_arr_st(min_value=np.float32(0.1), max_value=np.float32(1e5)))
    def test_comoving_distance(self, wrapper, vcosmo, z):
        """Test the wrapper's comoving_distance."""
        d = wrapper.comoving_distance(z)
        assert np.allclose(d, vcosmo.comoving_radial_distance(z, tol=0.0001))
        assert isinstance(d, np.ndarray)

    @settings(deadline=500)
    @given(z_arr_st(min_value=np.float32(0.1), max_value=np.float32(1e5)))
    def test_comoving_transverse_distance(self, wrapper, vcosmo, z):
        """Test the wrapper's comoving_transverse_distance."""
        d = wrapper.comoving_transverse_distance(z)
        assert np.allclose(d, vcosmo.angular_diameter_distance(z) * (1 + z))
        assert isinstance(d, np.ndarray)

    @settings(deadline=500)
    @given(z_arr_st())
    def test_comoving_volume_flat(self, wrapper, vcosmo, z):
        """Test the wrapper's comoving_volume."""
        v = wrapper.comoving_volume(z)
        assert isinstance(v, np.ndarray)
        if wrapper.Omega_k0 != 0:
            pytest.skip("comoving_volume_flat is only valid for flat cosmologies.")

        assert np.allclose(
            v, 4.0 / 3.0 * np.pi * vcosmo.comoving_radial_distance(z, tol=0.001) ** 3
        )

    @given(z_arr_st())
    def test_comoving_volume_positive(self, wrapper, cosmo, z):
        """Test the wrapper's comoving_volume."""
        v = wrapper.comoving_volume(z)
        assert isinstance(v, np.ndarray)
        if wrapper.Omega_k0 <= 0:
            pytest.skip(
                "comoving_volume_positive is only valid for positive cosmologies."
            )

        pytest.fail("TODO")

    @given(z_arr_st())
    def test_comoving_volume_negative(self, wrapper, cosmo, z):
        """Test the wrapper's comoving_volume."""
        v = wrapper.comoving_volume(z)
        assert isinstance(v, np.ndarray)
        if wrapper.Omega_k0 >= 0:
            pytest.skip(
                "comoving_volume_positive is only valid for negative cosmologies."
            )

        pytest.fail("TODO")

    @pytest.mark.skip("TODO")
    @given(z_arr_st())
    def test_differential_comoving_volume(self, wrapper, cosmo, z):
        """Test the wrapper's differential_comoving_volume."""
        v = wrapper.differential_comoving_volume(z)
        assert np.array_equal(v, cosmo.differential_comoving_volume(z))

        pytest.fail("TODO")

    @settings(deadline=500)
    @given(z_arr_st())
    def test_angular_diameter_distance(self, wrapper, vcosmo, z):
        """Test the wrapper's angular_diameter_distance."""
        d = wrapper.angular_diameter_distance(z)
        assert np.allclose(d, vcosmo.angular_diameter_distance(z))
        assert isinstance(d, np.ndarray)

    @settings(deadline=500)
    @given(z_arr_st(max_value=np.float32(1e5)))
    def test_luminosity_distance(self, wrapper, vcosmo, z):
        """Test the wrapper's luminosity_distance."""
        d = wrapper.luminosity_distance(z)
        assert np.allclose(d, vcosmo.luminosity_distance(z))
        assert isinstance(d, np.ndarray)
