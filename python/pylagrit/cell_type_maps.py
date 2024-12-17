from dataclasses import dataclass
from typing import cast, Any, List, Tuple

import numpy as np
import numpy.typing as npt

from pyvista import CellType, UnstructuredGrid


# Definition of number of nodes of element types
CELL_NODE_NUMBERS = {
    CellType.VERTEX: 1,
    CellType.LINE: 2,
    CellType.TRIANGLE: 3,
    CellType.QUAD: 4,
    CellType.TETRA: 4,
    CellType.PYRAMID: 5,
    CellType.WEDGE: 6,
    CellType.HEXAHEDRON: 8,
    None: 10,
}

# Node list of faces for element types
NODES_FOR_FACE = [
    None,  # point
    [[1], [2]],  # line
    [  # triangle
        [2, 3],
        [3, 1],
        [1, 2],
    ],
    [  # quadrilateral
        [1, 2],
        [2, 3],
        [3, 4],
        [4, 1],
    ],
    [  # tetrahedron
        [2, 3, 4],
        [1, 4, 3],
        [1, 2, 4],
        [1, 3, 2],
    ],
    [  # pyramid
        [1, 4, 3, 2],
        [1, 2, 5],
        [2, 3, 5],
        [3, 4, 5],
        [4, 1, 5],
    ],
    [  # prism
        [1, 3, 2],
        [4, 5, 6],
        [1, 2, 5, 4],
        [2, 3, 6, 5],
        [1, 4, 6, 3],
    ],
    [  # hexahedron
        [1, 4, 3, 2],
        [5, 6, 7, 8],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 4, 8, 7],
        [1, 5, 8, 4],
    ],
    None,  # hybrid
    None,  # polygon
]

# Cell type mapping to PFLOTRAN implicit unstructured grid
PF_CELL_TYPE_MAP = np.array(
    [
        None,  # vertex
        None,  # line
        None,  # triangle
        None,  # quad
        "T",  # tetra
        "P",  # pyramid
        "W",  # prism
        "H",  # hexahedron
        None,  # hybrid
        None,  # polygon
    ]
)

# Cell type mapping to VTK cell
VTK_CELL_TYPE_MAP = np.array(
    [
        CellType.VERTEX,
        CellType.LINE,
        CellType.TRIANGLE,
        CellType.QUAD,
        CellType.TETRA,
        CellType.PYRAMID,
        CellType.WEDGE,  # prism, VTK order = [0, 2, 1, 3, 5, 4]
        CellType.HEXAHEDRON,
        None,  # hybrid, no corresponding type
        None,  # polygon, no corresponding type
    ]
)


def lg_to_vtk_celltypes(itettyp: List[int]) -> npt.NDArray[Any]:
    return VTK_CELL_TYPE_MAP[np.array(itettyp, dtype=np.int64) - 1]


def vtk_to_lg_celltypes(cell_type: npt.NDArray[Any]) -> npt.NDArray[np.int64]:
    mapped_type = np.zeros_like(cell_type, dtype=np.int64)
    for i, ct in enumerate(cell_type):
        for j in range(8):
            if ct == VTK_CELL_TYPE_MAP[j]:
                mapped_type[i] = j + 1
                break

    return mapped_type


def lg_to_vtk_cells(
    itet: List[int], itettyp: List[int], itetoff: List[int]
) -> Tuple[npt.NDArray[np.int64], npt.NDArray[Any]]:
    el_nodes = np.array(itet) - 1
    el_offset = itetoff

    celltypes = lg_to_vtk_celltypes(itettyp)
    n_nodes = np.array([CELL_NODE_NUMBERS[t] for t in celltypes])
    cells = np.empty(len(n_nodes) + n_nodes.sum(), dtype=np.int64)

    el_idx = 0
    for offset, npts, typ in zip(el_offset, n_nodes, celltypes):
        cells[el_idx] = npts

        tmp = el_nodes[offset : offset + npts]
        if typ == CellType.WEDGE:  # prism type in lagrit
            tmp = tmp[[0, 2, 1, 3, 5, 4]]

        cells[el_idx + 1 : el_idx + 1 + npts] = tmp
        el_idx = el_idx + 1 + npts

    return cells, celltypes


@dataclass
class LgMeshData:
    xic: npt.NDArray[np.float64]
    yic: npt.NDArray[np.float64]
    zic: npt.NDArray[np.float64]
    itet: npt.NDArray[np.int64]
    itettyp: npt.NDArray[np.int64]
    itetoff: npt.NDArray[np.int64]

    @property
    def nnodes(self) -> int:
        return len(self.xic)

    @property
    def nelements(self) -> int:
        return len(self.itettyp)

    @property
    def nodes_per_element(self) -> int:
        vtk_cell_type = cast(CellType, VTK_CELL_TYPE_MAP[self.itettyp[0] - 1])
        return CELL_NODE_NUMBERS[vtk_cell_type]


def pv_ug_to_lg(ug: UnstructuredGrid) -> LgMeshData:
    xic = ug.points[:, 0]
    yic = ug.points[:, 1]
    zic = ug.points[:, 2]

    if len(ug.cells) == 0:
        return LgMeshData(xic, yic, zic, np.array([]), np.array([]), np.array([]))

    itettyp = vtk_to_lg_celltypes(ug.celltypes)

    itetoff = np.array([CELL_NODE_NUMBERS[t] for t in ug.celltypes]).cumsum()
    itetoff -= itetoff[0]

    itet = np.empty(
        itetoff[-1] + CELL_NODE_NUMBERS[CellType(ug.celltypes[-1])], dtype=np.int64
    )

    for i, offset in enumerate(itetoff):
        nnodes = CELL_NODE_NUMBERS[ug.celltypes[i]]

        tmp = ug.cells[offset + i + 1 : offset + i + 1 + nnodes] + 1
        if ug.celltypes[i] == CellType.WEDGE:  # prism type in lagrit
            tmp = tmp[[0, 2, 1, 3, 5, 4]]

        itet[offset : offset + nnodes] = tmp

    return LgMeshData(xic, yic, zic, itet, itettyp, itetoff)
