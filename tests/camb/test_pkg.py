"""Test the Cosmology API compat library."""


def test_imported():
    """This is a namespace package, so it should be importable."""
    import cosmology.compat.camb

    assert cosmology.compat.camb.__name__ == "cosmology.compat.camb"
