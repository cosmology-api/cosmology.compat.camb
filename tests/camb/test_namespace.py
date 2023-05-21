"""Test the Cosmology API compat library."""

import cosmology.compat.camb as namespace
from cosmology.api import CosmologyNamespace


def test_namespace_is_compliant():
    """Test :mod:`cosmology.compat.camb.constants`."""
    assert isinstance(namespace, CosmologyNamespace)
