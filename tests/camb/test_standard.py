"""Test the Cosmology API compat library."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from cosmology.api import (
    StandardCosmology,
    StandardCosmologyWrapper as StandardCosmologyWrapperAPI,
)

from .test_components import (
    BaryonComponent_Test,
    CurvatureComponent_Test,
    DarkEnergyComponent_Test,
    DarkMatterComponent_Test,
    MatterComponent_Test,
    NeutrinoComponent_Test,
    PhotonComponent_Test,
    TotalComponent_Test,
)
from .test_core import Test_CosmologyWrapper
from .test_distances import DistanceMeasures_Test
from .test_extras import CriticalDensity_Test, HubbleParameter_Test
from cosmology.compat.camb import StandardCosmologyWrapper

################################################################################
# TESTS
################################################################################


class Test_StandardCosmologyWrapper(
    TotalComponent_Test,
    CurvatureComponent_Test,
    MatterComponent_Test,
    BaryonComponent_Test,
    NeutrinoComponent_Test,
    DarkEnergyComponent_Test,
    DarkMatterComponent_Test,
    PhotonComponent_Test,
    CriticalDensity_Test,
    HubbleParameter_Test,
    DistanceMeasures_Test,
    Test_CosmologyWrapper,
):
    @pytest.fixture(scope="class")
    def wrapper(self, cosmo):
        return StandardCosmologyWrapper(cosmo)

    @pytest.fixture(scope="class")
    def vcosmo(self, cosmo):
        vc = SimpleNamespace()
        vc.comoving_radial_distance = np.vectorize(
            cosmo.comoving_radial_distance,
            excluded={"tol"},
        )
        vc.angular_diameter_distance = np.vectorize(cosmo.angular_diameter_distance)
        vc.luminosity_distance = np.vectorize(cosmo.luminosity_distance)
        vc.hubble_parameter = np.vectorize(cosmo.hubble_parameter)
        vc.get_Omega = np.vectorize(cosmo.get_Omega, excluded={"var"})
        return vc

    # =========================================================================
    # Tests

    def test_wrapper_is_compliant(self, wrapper):
        """Test that StandardCosmologyWrapper is a StandardCosmologyWrapper."""
        super().test_wrapper_is_compliant(wrapper)

        assert isinstance(wrapper, StandardCosmology)
        assert isinstance(wrapper, StandardCosmologyWrapperAPI)
