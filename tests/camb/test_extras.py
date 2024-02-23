"""Test the Cosmology API compat library."""

from __future__ import annotations

import numpy as np
from hypothesis import given

from cosmology.api import CriticalDensity, HubbleParameter

from .conftest import z_arr_st
from cosmology.compat.camb import constants

################################################################################
# TESTS
################################################################################


class CriticalDensity_Test:
    def test_wrapper_is_compliant(self, wrapper):
        """Test that AstropyCosmology is a BackgroundCosmologyWrapper."""
        if hasattr(super(), "test_wrapper_is_compliant"):
            super().test_wrapper_is_compliant(wrapper)

        assert isinstance(wrapper, CriticalDensity)

    def test_critical_density0(self, wrapper, params):
        """
        Test that the wrapper's critical_density0 is the same as
        critical_density0.
        """
        assert wrapper.critical_density0 == 3e6 * params.H0**2 / (
            8 * np.pi * constants.G
        )
        assert isinstance(wrapper.critical_density0, np.ndarray)

    @given(z_arr_st(max_value=1e9))
    def test_critical_density(self, wrapper, vcosmo, z):
        r"""Test that the wrapper's critical_density is critical_density."""
        rho = wrapper.critical_density(z)
        assert np.array_equal(
            rho, 3e6 * vcosmo.hubble_parameter(z) ** 2 / (8 * np.pi * constants.G)
        )
        assert isinstance(rho, np.ndarray)


class HubbleParameter_Test:
    def test_wrapper_is_compliant(self, wrapper):
        """Test that AstropyCosmology is a BackgroundCosmologyWrapper."""
        if hasattr(super(), "test_wrapper_is_compliant"):
            super().test_wrapper_is_compliant(wrapper)

        assert isinstance(wrapper, HubbleParameter)

    def test_H0(self, wrapper, params):
        """Test that the wrapper has the same H0 as the wrapped object."""
        assert wrapper.H0 == params.H0
        assert isinstance(wrapper.H0, np.ndarray)

    def test_hubble_distance(self, wrapper, params):
        """Test that the wrapper has the same hubble_distance as the wrapped object."""
        assert wrapper.hubble_distance == constants.c / params.H0
        assert isinstance(wrapper.hubble_distance, np.ndarray)

    def test_hubble_time(self, wrapper, params):
        """Test that the wrapper has the same hubble_time as the wrapped object."""
        assert wrapper.hubble_time == np.float64("978.5") / params.H0
        assert isinstance(wrapper.hubble_time, np.ndarray)

    @given(z_arr_st())
    def test_H(self, wrapper, vcosmo, z):
        """Test that the wrapper's H is the same as the wrapped object's."""
        H = wrapper.H(z)
        assert np.array_equal(H, vcosmo.hubble_parameter(z))
        assert isinstance(H, np.ndarray)

    @given(z_arr_st())
    def test_H_over_H0(self, wrapper, vcosmo, params, z):
        """Test that the wrapper's efunc is the same as the wrapped object's."""
        e = wrapper.H_over_H0(z)
        assert np.array_equal(e, vcosmo.hubble_parameter(z) / params.H0)
        assert isinstance(e, np.ndarray)
