"""The Cosmology API compatability wrapper for CAMB."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, cast

from cosmology.api import CosmologyNamespace, CosmologyWrapper
from numpy import floating
from numpy.typing import ArrayLike, NDArray

import camb

__all__: list[str] = []


@dataclass(frozen=True)
class CAMBCosmology(CosmologyWrapper[NDArray[floating[Any]], ArrayLike]):
    """The Cosmology API wrapper for :mod:`camb`."""

    cosmo: camb.CAMBdata
    name: str | None = None

    def __post_init__(self) -> None:
        """Run-time post-processing.

        Note that if this module is c-compiled (e.g. with :mod:`mypyc`) that
        the type of ``self.cosmo`` must be ``CAMBdata`` at object creation
        and cannot be later processed here.
        """
        if not isinstance(self.cosmo, (camb.CAMBdata, camb.CAMBParams)):
            msg = (
                "cosmo must be a CAMBdata or CAMBParams instance, "
                f"not {type(self.cosmo)}"
            )
            raise TypeError(msg)
        elif isinstance(self.cosmo, camb.CAMBParams):
            cosmo = camb.get_background(self.cosmo)
            params = self.cosmo

            # Need to fix the type of cosmo to CAMBdata.
            object.__setattr__(self, "cosmo", cosmo)
        else:
            params = self.cosmo.Params

        self.params: camb.CAMBparams
        object.__setattr__(self, "params", params)

    @property
    def __cosmology_namespace__(self, /) -> CosmologyNamespace:
        """Returns a `CosmologyNamespace` with the cosmology API functions."""
        import cosmology.compat.camb as namespace

        return cast(CosmologyNamespace, namespace)
