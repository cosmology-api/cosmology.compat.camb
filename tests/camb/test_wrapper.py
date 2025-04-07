import astropy.cosmology
import camb
import numpy as np
import pytest

import cosmology.compat.camb


@pytest.fixture(scope="session")
def rng():
    return np.random.default_rng(42)


@pytest.fixture(scope="session")
def z(rng):
    return np.sort(rng.uniform(0.0, 5.0, 100))


@pytest.fixture(scope="session")
def spaced_z(rng):
    r"""Return $z \in [0, 500]$ where $z_{i+1} - z_{i} >= 1$."""
    # Settings
    low, high = 0, 500
    min_spacing = 1
    n = 100

    # Generate random z values
    total_gap = (high - low) - min_spacing * (n - 1)
    assert total_gap > 0, "Total gap must be positive"
    # sample random "extra" gaps to add between each pair
    extra_gaps = rng.dirichlet(np.ones(n - 1)) * total_gap
    gaps = min_spacing + extra_gaps
    # cumulative sum to get positions
    return low + np.concatenate([[0], np.cumsum(gaps)])


@pytest.fixture(scope="session")
def compare():
    return astropy.cosmology.LambdaCDM(
        H0=70.0,
        Om0=0.3,
        Ob0=0.05,
        Ode0=0.8,
        Tcmb0=2.73,
        Neff=0.0,
        m_nu=0.0,
    )


@pytest.fixture(scope="session")
def cosmo(compare):
    pars = camb.set_params(
        H0=compare.H0.value,
        omch2=(compare.Om0 - compare.Ob0) * compare.h**2,
        ombh2=compare.Ob0 * compare.h**2,
        omnuh2=0.0,
        omk=compare.Ok0,
        TCMB=compare.Tcmb0.value,
        num_nu_massless=0.0,
        num_nu_massive=0,
        mnu=0.0,
        nnu=0.0,
    )
    results = camb.get_background(pars)
    return cosmology.compat.camb.Cosmology(pars, results)


def test_h(cosmo, compare):
    assert cosmo.h == compare.h


def test_H0(cosmo, compare):
    assert compare.H0.value == cosmo.H0


def test_Omega_m0(cosmo, compare):
    np.testing.assert_allclose(
        cosmo.Omega_m0,
        compare.Om0,
        rtol=1e-10,
    )


def test_Omega_de0(cosmo, compare):
    np.testing.assert_allclose(
        cosmo.Omega_de0,
        compare.Ode0,
        rtol=1e-10,
    )


def test_Omega_k0(cosmo, compare):
    assert cosmo.Omega_k0 == compare.Ok0


def test_hubble_distance(cosmo, compare):
    np.testing.assert_allclose(
        cosmo.hubble_distance,
        compare.hubble_distance.value,
        rtol=1e-10,
    )


def test_H(z, cosmo, compare):
    np.testing.assert_allclose(
        cosmo.H(z),
        compare.H(z).value,
        rtol=1e-12,
    )


def test_Omega_m(z, cosmo, compare):
    np.testing.assert_allclose(
        cosmo.Omega_m(z),
        compare.Om(z) + compare.Onu(z),
        rtol=1e-3,
    )


def test_Omega_de(z, cosmo, compare):
    np.testing.assert_allclose(
        cosmo.Omega_de(z),
        compare.Ode(z),
        rtol=1e-3,
    )


def test_Omega_k(z, cosmo, compare):
    np.testing.assert_allclose(
        cosmo.Omega_k(z),
        compare.Ok(z),
        rtol=1e-3,
    )


def test_comoving_distance(z, cosmo, compare):
    np.testing.assert_allclose(
        cosmo.comoving_distance(z),
        compare.comoving_distance(z).value,
        rtol=1e-3,
    )

    z1, z2 = z[:-1], z[1:]

    np.testing.assert_allclose(
        cosmo.comoving_distance(z1, z2),
        compare.comoving_distance(z2).value - compare.comoving_distance(z1).value,
        rtol=1e-3,
    )


def test_inv_comoving_distance(z, cosmo, compare):
    x = compare.comoving_distance(z).value
    np.testing.assert_allclose(
        cosmo.inv_comoving_distance(x),
        z,
        rtol=1e-3,
    )


def test_angular_diameter_distance(z, cosmo, compare):
    np.testing.assert_allclose(
        cosmo.angular_diameter_distance(z),
        compare.angular_diameter_distance(z).value,
        rtol=1e-3,
    )

    z1, z2 = z[:-1], z[1:]

    np.testing.assert_allclose(
        cosmo.angular_diameter_distance(z1, z2),
        compare.angular_diameter_distance_z1z2(z1, z2).value,
        rtol=1e-3,
    )


def test_H_over_H0(z, cosmo, compare):
    np.testing.assert_allclose(
        cosmo.H_over_H0(z),
        compare.efunc(z),
        rtol=1e-12,
    )


def test_transverse_comoving_distance(spaced_z, cosmo, compare):
    np.testing.assert_allclose(
        cosmo.transverse_comoving_distance(spaced_z),
        compare.comoving_transverse_distance(spaced_z).value,
        rtol=1e-4,
    )

    z1, z2 = spaced_z[:-1], spaced_z[1:]

    np.testing.assert_allclose(
        cosmo.transverse_comoving_distance(z1, z2),
        compare.comoving_transverse_distance(z2).value
        - compare.comoving_transverse_distance(z1).value,
        rtol=1,
    )
