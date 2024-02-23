"""Test the Cosmology API compat library."""

from cosmology.api import CosmologyNamespace

import cosmology.compat.camb as namespace


def test_namespace_is_compliant():
    """Test :mod:`cosmology.compat.camb.constants`."""
    assert isinstance(namespace, CosmologyNamespace)
