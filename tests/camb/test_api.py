from typing import Any, Protocol, runtime_checkable

import camb
from numpy.typing import NDArray

import cosmology.api
import cosmology.compat.camb

Redshift = float | NDArray[Any]
Array = NDArray[Any]


@runtime_checkable
class Cosmology(
    cosmology.api.HasH0,
    cosmology.api.HasOmegaM0,
    cosmology.api.HasOmegaDE0,
    cosmology.api.HasOmegaK0,
    cosmology.api.HasHubbleDistance,
    cosmology.api.HasOmegaM,
    cosmology.api.HasOmegaDE,
    cosmology.api.HasOmegaK,
    cosmology.api.HasComovingDistance,
    cosmology.api.HasAngularDiameterDistance,
    Protocol,
): ...


def test_api():
    pars = camb.set_params(H0=70.0)
    results = camb.get_background(pars)

    cosmo = cosmology.compat.camb.Cosmology(pars, results)

    assert isinstance(cosmo, Cosmology)
