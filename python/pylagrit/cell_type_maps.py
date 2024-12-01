import numpy as np

from pyvista import CellType


# Definition of number of nodes for element type
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
