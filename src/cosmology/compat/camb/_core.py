"""The Cosmology API compatibility wrapper for CAMB."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Union, cast

import camb
from numpy import floating
from numpy.typing import NDArray

from cosmology.api import CosmologyNamespace
from cosmology.api.compat import CosmologyWrapper as CosmologyWrapperAPI

if TYPE_CHECKING:
    from typing_extensions import TypeAlias


NDFloating: TypeAlias = NDArray[floating[Any]]
InputT: TypeAlias = Union[NDFloating, float]


@dataclass(frozen=True)
class CosmologyWrapper(CosmologyWrapperAPI[NDFloating, InputT]):
    """The Cosmology API wrapper for :mod:`camb`."""

    cosmo: camb.CAMBdata
    name: str | None = None

    def __post_init__(self) -> None:
        """Run-time post-processing.

        Note that if this module is c-compiled (e.g. with :mod:`mypyc`) that
        the type of ``self.cosmo`` must be ``CAMBdata`` at object creation
        and cannot be later processed here.
        """
        if not isinstance(self.cosmo, (camb.CAMBdata, camb.CAMBparams)):
            msg = (
                "cosmo must be a CAMBdata or CAMBParams instance, "
                f"not {type(self.cosmo)}"
            )
            raise TypeError(msg)
        elif isinstance(self.cosmo, camb.CAMBparams):
            cosmo = camb.get_background(self.cosmo)
            params = self.cosmo

            # Need to fix the type of cosmo to CAMBdata.
            object.__setattr__(self, "cosmo", cosmo)
        else:
            params = self.cosmo.Params

        self._params: camb.CAMBparams
        object.__setattr__(self, "_params", params)

        self._cosmo_fn: dict[str, Any]
        object.__setattr__(self, "_cosmo_fn", {})

    @property
    def __cosmology_namespace__(self) -> CosmologyNamespace:
        """Returns an object that has all the cosmology API functions on it.

        Returns
        -------
        `CosmologyNamespace`
            An object representing the CAMB cosmology API namespace.

        """
        import cosmology.compat.camb as namespace

        return cast(CosmologyNamespace, namespace)
