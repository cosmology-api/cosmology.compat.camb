"""Test the Cosmology API compat library."""

from types import SimpleNamespace

import camb
import pytest

from cosmology.api import (
    Cosmology as CosmologyAPI,
    CosmologyWrapper as CosmologyWrapperAPI,
)

from cosmology.compat.camb._core import CosmologyWrapper

################################################################################
# TESTS
################################################################################


class Test_CosmologyWrapper:
    @pytest.fixture(scope="class")
    def cosmo(self) -> camb.CAMBdata:
        # EXAMPLE FROM https://camb.readthedocs.io/en/latest/CAMBdemo.html

        # Set up a new set of parameters for CAMB
        pars = camb.CAMBparams()
        # This function sets up CosmoMC-like settings, with one massive neutrino
        # and helium set using BBN consistency
        pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
        pars.InitPower.set_params(As=2e-9, ns=0.965, r=0)
        pars.set_for_lmax(2500, lens_potential_accuracy=0)

        return camb.get_results(pars)

    @pytest.fixture(scope="class")
    def vcosmo(self, cosmo):
        return SimpleNamespace()

    @pytest.fixture(scope="class")
    def params(self, cosmo):
        return cosmo.Params

    @pytest.fixture(scope="class")
    def wrapper(self, cosmo):
        return CosmologyWrapper(cosmo)

    # =========================================================================
    # Tests

    def test_wrapper_is_compliant(self, wrapper):
        """Test that CosmologyWrapper is a CosmologyWrapper."""
        if hasattr(super(), "test_wrapper_is_compliant"):
            super().test_wrapper_is_compliant(wrapper)

        assert isinstance(wrapper, CosmologyAPI)
        assert isinstance(wrapper, CosmologyWrapperAPI)

    def test_getattr(self, wrapper, cosmo):
        """Test that the wrapper can access the attributes of the wrapped object."""
        assert wrapper.grhocrit == cosmo.grhocrit
