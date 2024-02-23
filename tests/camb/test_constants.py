"""Test the Cosmology API compat library."""

from cosmology.api import CosmologyConstantsNamespace

from cosmology.compat.camb import constants


def test_namespace_is_compliant():
    """Test :mod:`cosmology.compat.camb.constants`."""
    assert isinstance(constants, CosmologyConstantsNamespace)
