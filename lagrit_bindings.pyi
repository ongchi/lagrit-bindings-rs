import typing

import numpy as np
import numpy.typing as npt

IntType = typing.TypeVar("IntType", bound=np.integer)
IntArray = typing.Annotated[npt.NDArray[IntType], typing.Literal["N"]]

FloatType = typing.TypeVar("FloatType", bound=np.float64)
FloatArray = typing.Annotated[npt.NDArray[FloatType], typing.Literal["N"]]

class MeshObject:
    def __init__(self, name: str): ...
    def name(self) -> str: ...
    def status(self): ...
    def mesh_type(self) -> typing.Tuple[str, int]: ...
    def attr(
        self, attr: str
    ) -> (
        int | typing.List[int] | float | typing.List[float] | str | typing.List[str]
    ): ...
    def set_attr(self, attr: str, value: IntArray | FloatArray): ...
    def attr_list(self) -> typing.List[typing.List[str]]: ...
    def pset_names(self) -> typing.List[str]: ...
    def eltset_names(self) -> typing.List[str]: ...

class LaGriT:
    def __init__(
        self,
        mode: str,
        log_file: typing.Optional[str],
        batch_file: typing.Optional[str],
        workdir: typing.Optional[str],
    ): ...
    def sendcmd(self, cmd: str): ...
    def cmdmsg(self) -> str: ...
    def mo_names(self) -> typing.List[str]: ...
    def pset_names(self) -> typing.List[str]: ...
    def eltset_names(self) -> typing.List[str]: ...
    def cmo(self) -> MeshObject: ...
    def close(self): ...
