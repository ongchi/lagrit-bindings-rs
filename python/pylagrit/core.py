from collections import defaultdict
from itertools import product
import xml.etree.ElementTree as ET

from pathlib import Path
from typing import (
    Annotated,
    Dict,
    Iterable,
    List,
    Literal,
    Optional,
    OrderedDict,
    Tuple,
    TypeVar,
    cast,
)
from xml.dom import minidom

import numpy as np
import numpy.typing as npt
import pyvista as pv
from pyvista import CellType

from lagrit_bindings import LaGriT, MeshObject
from pylagrit.cell_type_maps import (
    CELL_NODE_NUMBERS,
    NODES_FOR_FACE,
    PF_CELL_TYPE_MAP,
    VTK_CELL_TYPE_MAP,
)

DType = TypeVar("DType", bound=np.generic)
Points = Annotated[npt.NDArray[DType], Literal["N", 3]]


def new_name(base, names):
    i = 1
    name = base + str(i)
    while name in names:
        i += 1
        name = base + str(i)
    return name


# TODO: LaGriT does no alwasy returns an error when something goes wrong.
#       Find error messages in output file or stdout to capture all errors.
class PyLaGriT:
    """
    Python lagrit class

    Parameters
    ----------
    mode : str, optional
        tty output mode, by default "noisy"

    log_file : str, optional
        Path to log file, by default "lagrit.log"

    out_file : str, optional
        Path to out file, by default "lagrit.out"

    """

    def __init__(
        self,
        mode: Literal["slient", "noisy"] = "noisy",
        log_file: str = "lagrit.log",
        out_file: str = "lagrit.out",
    ):
        self.core = LaGriT(mode, log_file, out_file)

    def sendcmd(self, cmd: str):
        self.core.sendcmd(cmd)

    def close(self):
        self.core.close()

    @property
    def before(self) -> str:
        return self.core.cmdmsg()

    @property
    def mo_list(self) -> List["MO"]:
        return [MO(k, self) for k in self.core.mo_names()]

    @property
    def cmo(self) -> "MO":
        cmo = self.core.cmo()
        return MO(cmo.name(), self)

    def read_mo(
        self,
        filename: str,
        filetype: Optional[str] = None,
        name: Optional[str] = None,
        binary=False,
    ) -> Optional["MO" | List["MO"]]:
        """
        Read in mesh

        Parameters
        ----------
        filename : str
            Name of mesh file to read in

        filetype : str, optional
            Type of file, automatically detected if not specified, by default None

        name : str, optional
            Internal Lagrit name of new mesh object, automatically created if None, by default None

        binary : bool, optional
            Indicates that file is binary if True, ascii if False, by default False

        Returns
        -------
        MO

        Examples
        --------
        Example 1

        >>> import pylagrit
        >>> import numpy as np
        >>> lg = pylagrit.PyLaGriT()

        Create a mesh object and dump it to a gmv file 'test.gmv'.
        >>> mo = lg.create(name="test")
        >>> mo.createpts_brick("xyz", (5, 5, 5), (0, 0, 0), (5, 5, 5))
        >>> mo.dump("test.gmv")
        >>> mo.dump("test.avs")

        >>> mo1 = lg.read_mo("test.gmv")
        >>> mo2 = lg.read_mo("test.avs")

        >>> lg.close()
        >>> lg = pylagrit.PyLaGriT()

        Example 2 - Reading in LaGriT binary file

        Create list with mesh object as first element
        >>> dxyz = (0.25, 0.25, 0.25)
        >>> mins = (0.0, 0.0, 0.0)
        >>> maxs = (1.0, 1.0, 1.0)
        >>> mos = [
        ...     lg.createpts_dxyz(dxyz, mins, maxs, "tet", connect=True, name="testmo")
        ... ]

        Create three new mesh objects, each one directly above the other
        >>> for i in range(3):
        ...     mos.append(mos[-1].copy())
        ...     mins = (mos[-1].attr("xmin"), mos[-1].attr("ymin"), mos[-1].attr("zmin"))
        ...     mos[-1].trans(mins, mins + np.array([0.0, 0.0, 1.0]))
        >>> lg.dump("lagrit_binary.lg")

        >>> lg.close()
        >>> lg = pylagrit.PyLaGriT()

        >>> ms_read = lg.read_mo("lagrit_binary.lg")
        >>> print([mo.name for mo in ms_read])
        ['testmo', 'mo1', 'mo2', 'mo3']

        >>> lg.close()

        """
        old_mo_list = self.core.mo_names()

        cmd = ["read"]

        # If filetype is lagrit, name is irrelevant
        islg = filetype == "lagrit" or filename.split(".")[-1].lower() in [
            "lg",
            "lagrit",
        ]

        if islg:
            cmd.append("lagrit")

        cmd.append(filename)

        if filetype is not None:
            cmd.append(filetype)

        if not islg and name is None:
            name = new_name("mo", self.core.mo_names())
            cmd.append(name)

        if binary:
            cmd.append("binary")

        self.sendcmd("/".join(cmd))

        if islg:
            cur_mo_list = self.core.mo_names()
            for mo in old_mo_list:
                cur_mo_list.remove(mo)
            return [MO(mo, self) for mo in cur_mo_list]
        else:
            return MO(cast(str, name), self)

    def read_fehm(
        self, filename: str, avs_filename="temp.inp", elem_type: Optional[str] = None
    ):
        with open(filename) as fh:
            ln = fh.readline()
            nn = int(fh.readline().strip())
            while "elem" not in ln:
                ln = fh.readline()
            vs = fh.readline().strip().split()
        elem_int = int(vs[0])
        ne = int(vs[1])
        crds = np.genfromtxt(filename, skip_header=2, max_rows=nn)
        conns = np.genfromtxt(filename, skip_header=2 + nn + 3, max_rows=ne)
        with open(avs_filename, "w") as fh:
            fh.write("    %d    %d    0    0    0\n" % (nn, ne))
            np.savetxt(fh, crds, fmt="%d %f %f %f")
            if elem_type is None:
                if elem_int == 8:
                    elem_type = "hex"
                elif elem_int == 3:
                    elem_type = "tri"
                elif elem_int == 4:
                    if (
                        np.all(np.diff(crds[:, 1])) == 0
                        or np.all(np.diff(crds[:, 2])) == 0
                        or np.all(np.diff(crds[:, 3])) == 0
                    ):
                        elem_type = "qua"
                    else:
                        elem_type = "tet"
            for conn in conns:
                fh.write("%d 1 %s" % (conn[0], elem_type))
                for i in range(elem_int):
                    fh.write(" %d" % conn[i + 1])
                fh.write("\n")
        return self.read_mo(avs_filename)

    def read_sheetij(
        self,
        name: str,
        filename: str,
        NXY: Tuple[int, int],
        minXY: Tuple[float, float],
        DXY: Tuple[float, float],
        connect=True,
        file_type: Literal["ascii", "binary"] = "ascii",
        flip: Optional[Literal["x", "y", "xy"]] = None,
        skip_lines=0,
        data_type: Literal["float", "double"] = "float",
    ) -> "MO":
        """
        Creates a quad mesh from an elevation file. Note the input file is read as Z(i,j) into the cmo attribute 'zic'

        Parameters
        ----------
        name : str
            Name of mesh object

        filename : str
            Elevation filename

        NXY : Tuple[int, int]
            [nx, ny] - [columns in x-direction, rows in y-direction]

        minXY : Tuple[float, float]
            [minX, minY] - location of lower left corner

        DXY : Tuple[float, float]
            [Dx, Dy] - cell size in x and y directions

        connect : bool, optional
            True will create a quad grid, otherwise keeps data as points, by default True

        file_type : str, optional
            May be either ascii or binary, by default "ascii"

        flip : str, optional
            May be 'x', 'y' to reflect across those axes, or 'none' to keep static, by default "none"

        skip_lines : int, optional
            Skip n number of header lines, by default 0

        data_type : str, optional
            Read in elevation data as either float or double, by default "float"

        Returns
        -------
        MO

        Examples
        --------
        Building a surface mesh from Modflow elevation file

        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()
        >>>
        >>> # Elevation files are typically headerless unwrapped vectors
        >>> # Define parameters to pack these elements into a matrix
        >>> ncols = 276
        >>> nrows = 313
        >>> DXY = (100, 100)

        >>> elev_surface = lg.read_sheetij(
        ...     "surfacemesh", "example.mod", [ncols, nrows], [0, 0], DXY, flip="y"
        ... )

        >>> lg.close()

        """

        connect_str = "connect" if connect else "points"
        skip_str = "skip %d" % skip_lines

        if flip is None:
            flip_str = ""
        elif flip == "x":
            flip_str = "xflip"
        elif flip == "y":
            flip_str = "yflip"
        elif flip == "xy":
            flip_str = "xflip,yflip"

        # Create new mesh object with given name
        self.sendcmd(f"cmo/create/{name}")
        self.sendcmd(f"cmo/select/{name}")

        # Read in elevation file and append to mesh
        cmd = [
            "read",
            "sheetij",
            filename,
            ",".join([str(v) for v in NXY]),
            ",".join([str(v) for v in minXY]),
            ",".join([str(v) for v in DXY]),
            skip_str,
            flip_str,
            connect_str,
            file_type,
            data_type,
        ]
        self.sendcmd("/".join(cmd))

        return MO(name, self)

    def boundary_components(
        self,
        style: Literal["node", "element"] = "node",
        material_id_number: Optional[int] = None,
        reset: Optional[bool] = None,
    ):
        """
        Calculates the number of connected components of a mesh for diagnostic purposes.

        Parameters
        ----------
        style : str, optional
            May be either 'node' or 'element', by default 'node'

        material_id_number : int, optional
            Only examines nodes with imt = mat. id number, by default None

        reset : bool, optional
            May be either True, False, or None, by default None

        """

        cmd = ["boundary_components", style]

        if material_id_number:
            cmd.append(str(material_id_number))
        if reset is not None:
            if reset:
                cmd.append("reset")
            elif not reset:
                cmd.append("noreset")

        self.sendcmd("/".join(cmd))

    def addmesh(
        self,
        mo1: "MO | str",
        mo2: "MO | str",
        style="add",
        name: Optional[str] = None,
        *args,
    ):
        if name is None:
            name = new_name("mo", self.core.mo_names())
        cmd = ["addmesh", style, name, str(mo1), str(mo2)]
        for a in args:
            if isinstance(a, str):
                cmd.append(a)
            elif isinstance(a, list):
                cmd.extend([str(v) for v in a])
        self.sendcmd("/".join(cmd))
        return MO(name, self)

    def addmesh_add(
        self,
        mo1: "MO | str",
        mo2: "MO | str",
        name: Optional[str] = None,
        refine_factor=-1,
        refine_style="edge",
    ):
        return self.addmesh(mo1, mo2, "add", name, refine_factor, refine_style)

    def addmesh_amr(self, mo1: "MO | str", mo2: "MO | str", name: Optional[str] = None):
        return self.addmesh(mo1, mo2, style="amr", name=name)

    def addmesh_append(
        self, mo1: "MO | str", mo2: "MO | str", name: Optional[str] = None
    ):
        return self.addmesh(mo1, mo2, style="append", name=name)

    def addmesh_delete(
        self, mo1: "MO | str", mo2: "MO | str", name: Optional[str] = None
    ):
        return self.addmesh(mo1, mo2, style="delete", name=name)

    def addmesh_glue(
        self, mo1: "MO | str", mo2: "MO | str", name: Optional[str] = None
    ):
        return self.addmesh(mo1, mo2, style="glue", name=name)

    def addmesh_intersect(
        self,
        pset: "PSet | str",
        mo1: "MO | str",
        mo2: "MO | str",
        name: Optional[str] = None,
    ):
        if name is None:
            name = new_name("mo", self.core.mo_names())
        cmd = "/".join(["addmesh", "intersect", name, str(pset), str(mo1), str(mo2)])
        self.sendcmd(cmd)
        return PSet(name, self)

    def addmesh_merge(
        self, mo1: "MO | str", mo2: "MO | str", name: Optional[str] = None
    ):
        return self.addmesh(mo1, mo2, style="merge", name=name)

    def addmesh_pyramid(
        self, mo1: "MO | str", mo2: "MO | str", name: Optional[str] = None
    ):
        return self.addmesh(mo1, mo2, style="pyramid", name=name)

    def addmesh_excavate(
        self,
        mo1: "MO | str",
        mo2: "MO | str",
        name: Optional[str] = None,
        bfs=False,
        connect=False,
    ):
        if bfs:
            bfsstr = "bfs"
        else:
            bfsstr = " "
        if connect:
            connectstr = "connect"
        else:
            connectstr = " "
        return self.addmesh(mo1, mo2, "excavate", name, bfsstr, connectstr)

    def extract_surfmesh(
        self,
        name: Optional[str] = None,
        cmo_in: "Optional[MO | str]" = None,
        stride=(1, 0, 0),
        reorder=True,
        resetpts_itp=True,
        external=False,
        append=None,
    ) -> "MO":
        if name is None:
            name = new_name("mo", self.core.mo_names())

        stride = [str(v) for v in stride]
        cmd = ["extract/surfmesh", ",".join(stride), name]

        if cmo_in is not None:
            if isinstance(cmo_in, str):
                cmo = MO(cmo_in, self)
            else:
                cmo = cmo_in

            if resetpts_itp:
                cmo.resetpts_itp()

            if reorder:
                cmo.sendcmd("createpts/median")
                self.sendcmd(
                    "/".join(
                        [
                            "sort",
                            str(cmo_in),
                            "index/ascending/ikey/itetclr zmed ymed xmed",
                        ]
                    )
                )
                self.sendcmd("/".join(["reorder", str(cmo_in), "ikey"]))
                self.sendcmd("/".join(["cmo/DELATT", str(cmo_in), "xmed"]))
                self.sendcmd("/".join(["cmo/DELATT", str(cmo_in), "ymed"]))
                self.sendcmd("/".join(["cmo/DELATT", str(cmo_in), "zmed"]))
                self.sendcmd("/".join(["cmo/DELATT", str(cmo_in), "ikey"]))

            cmd.append(str(cmo_in))

        if external:
            cmd.append("external")

        if append:
            cmd.append(append)

        self.sendcmd("/".join(cmd))
        return MO(name, self)

    def read_script(self, fname: str):
        """
        Read a LaGriT Script

        Given a script name, executes the script in LaGriT.

        Parameters
        ----------
        fname : str
            The name or path to the lagrit script.

        """

        f = open(fname)
        commands = f.readlines()
        for c in commands:
            # Remove newlines and spaces
            c = "".join(c.split())
            if len(c) != 0 and "finish" not in c:
                self.sendcmd(c)

    def read_att(
        self,
        fname: str,
        attributes: Iterable[str],
        mesh: Optional["MO"] = None,
        operation="add",
    ):
        """
        Reads data from a file into an attribute.
        """

        if mesh is None:
            mesh = self.create()

        if isinstance(operation, (list, tuple)):
            operation = ",".join(list(map(str, operation)))

        cmd = "/".join(
            ["cmo", "readatt", mesh.name, ",".join(attributes), operation, fname]
        )
        self.sendcmd(cmd)

        return mesh

    def define(self, **kwargs):
        """
        Pass in a variable number of arguments to be defined in
        LaGriT's internal global scope.

        Note that it is generally considered bad practice in PyLaGriT
        to rely on LaGriT's variable system for parameters; however,
        there are use-cases where it is necessary: i.e., macro scripts.

        Examples
        --------

        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()

        >>> lg.define(MO_PTS="mo1",OUTFILE="mesh.inp",PERTURB32=1.3244)

        >>> lg.close()

        Equivalent to LaGriT script:

        define / MO_PTS / mo1
        define / OUTFILE / mesh.inp
        define / PERTURB32 / 1.3244

        """

        for key, value in kwargs.items():
            self.sendcmd(f"define / {key} / {value}")

    def convert(self, filename: str, to_type: str):
        """
        Convert File

        For each file of the pattern, creates a new file in the to_format format.
        The new files will be inside the directory that the LaGriT object was
        instantiated. The name of each file will be the same as the original
        file with the extension changed to to_type.

        Supports conversion from avs, and gmv files.
        Supports conversion to avs, exo, and gmv files.

        Parameters
        ----------
        filename : str
            Name of the file to be converted.

        to_type : str
            New format to convert files.

        Examples
        --------
        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()

        Create a mesh object and dump it to a gmv file 'test.gmv'.
        >>> mo = lg.create(name="test")
        >>> mo.createpts_brick(
        ...     "xyz",
        ...     (5, 5, 5),
        ...     (0, 0, 0),
        ...     (5, 5, 5),
        ... )
        >>> mo.dump("gmv", "test.gmv")

        Convert test.gmv to exoduce and contour files.
        >>> lg.convert("test.gmv", "exo")
        >>> lg.convert("test.gmv", "avs")

        >>> lg.close()

        """
        # Make sure I support the new filetype.
        if to_type not in ["avs", "gmv", "exo"]:
            raise ValueError(f"Conversion to {to_type} not supported.")

        fname = Path(filename)
        # Check that I support the old filetype.
        if fname.suffix not in [".avs", ".gmv"]:
            raise ValueError(f"Conversion from {fname.suffix} not supported.")

        cmo: MO = self.read_mo(str(filename))  # type: ignore
        cmo.dump(f"{fname.stem}.{to_type}")
        cmo.delete()

    def merge(self, mesh_objs: List["MO"], name: Optional[str] = None) -> "MO":
        """
        Merge Mesh Objects

        Merges two or more mesh objects together and returns the combined mesh
        object.

        Parameters
        ----------
        mesh_objs : List[MO]
            List of mesh objects to merge.

        Returns
        -------
        MO

        Examples
        --------
        >>> import pylagrit
        >>> import numpy as np
        >>> lg = pylagrit.PyLaGriT()

        Create list with mesh object as first element
        >>> dxyz = (0.25, 0.25, 0.25)
        >>> mins = (0.0, 0.0, 0.0)
        >>> maxs = (1.0, 1.0, 1.0)
        >>> mos = [lg.createpts_dxyz(dxyz, mins, maxs, "tet", connect=True)]

        Create three new mesh objects, each one directly above the other
        >>> for i in range(3):
        ...     mos.append(mos[-1].copy())
        ...     mins = (mos[-1].attr("xmin"), mos[-1].attr("ymin"), mos[-1].attr("zmin"))
        ...     mos[-1].trans(mins, mins + np.array([0.0, 0.0, 1.0]))

        Merge list of mesh objects and clean up
        >>> mo_merge = lg.merge(mos)
        >>> for mo in mos:
        ...     mo.delete()
        >>> mo_merge.rmpoint_compress(filter_bool=True, resetpts_itp=True)

        >>> lg.close()

        """
        if name is None:
            name = new_name("mo", self.core.mo_names())
        if len(mesh_objs) > 1:
            for mo in mesh_objs:
                cmd = "/".join(["addmesh", "merge", name, name, mo.name])
                self.sendcmd(cmd)
        else:
            raise ValueError("Must provide at least two objects to merge.")
        return MO(name, self)

    def create(
        self,
        elem_type: Literal[
            "tet", "hex", "pri", "pyr", "tri", "qua", "hyb", "lin", "triplane"
        ] = "tet",
        name: Optional[str] = None,
        npoints=0,
        nelements=0,
    ) -> "MO":
        """
        Create a Mesh Object

        Creates a mesh object in lagrit and an MO in the LaGriT object. Returns
        the mesh object.

        Parameters
        ----------
        elem_type : str, optional
            The type of mesh object to create, by default "tet"
            Options include:
                tet: Tetrahedron
                hex: Hexagon
                pri: Prism
                pyr: Pyramid
                tri: Triangle
                qua: Quadrilateral
                hyb: Hybrid
                lin: Line
                triploane: Triplane

        name : str, optional
            Name to be given to the mesh object, by default None

        npoints : int, optional
            The number of points in the mesh object, by default 0

        nelements : int, optional
            The number of elements in the mesh object, by default 0

        Returns
        -------
        MO

        """

        if name is None:
            name = new_name("mo", self.core.mo_names())

        self.sendcmd("cmo/create/%s/%i/%i/%s" % (name, npoints, nelements, elem_type))
        return MO(name, self)

    def create_line(
        self,
        npoints=0,
        mins: Optional[Tuple[float, float, float]] = None,
        maxs: Optional[Tuple[float, float, float]] = None,
        rz_switch=(1, 1, 1),
        name: Optional[str] = None,
    ) -> "MO":
        """Create a line mesh object."""
        mo_new = self.create(elem_type="lin", name=name, npoints=npoints)
        if mins is not None and maxs is not None:
            mo_new.createpts_line(npoints, mins, maxs, rz_switch)
        return mo_new

    def copy(self, mo: "MO | str", name: Optional[str] = None) -> "MO":
        """
        Copy Mesh Object

        Copies a mesh object, mo, and returns the MO object.

        Parameters
        ----------
        mo : MO or str
            Mesh object to copy.

        name : str, optional
            Name of the new mesh object, by default None

        Returns
        -------
        MO
        """

        # Check if name was specified, if not just generate one.
        if name is None:
            name = new_name("mo", self.core.mo_names())

        # Create the MO in lagrit and the PyLaGriT object.
        self.sendcmd(f"cmo/copy/{name}/{str(mo)}")
        return MO(name, self)

    def dump(
        self,
        filename: str,
        mos: "Optional[Iterable[MO | str]]" = None,
        filetype: Literal["binary", "ascii"] = "binary",
    ):
        """
        Dump lagrit binary file

        Parameters
        ----------
        filename : str
            Name of lagrit binary file to create

        mos : Iterable[MO]
            List of mesh objects to include

        filetype : str, optional
            Filetype to dump, 'binary' or 'ascii', by default "binary"

        """
        cmd = ["dump", "lagrit", filename]
        if mos is None:
            cmd.append("-all-")
        else:
            cmd.extend([str(mo) for mo in mos])
        if filetype == "ascii":
            cmd.append("ascii")
        self.sendcmd("/".join(cmd))

    def tri_mo_from_polyline(
        self,
        coords: List[Tuple[float, float, float] | Tuple[float, float]],
        filename="polyline.inp",
        name: Optional[str] = None,
    ) -> "MO":
        """
        Create polygon tri mesh object from points
        Points are expected to be defined clockwise by default

        Parameters
        ----------
        coords : List[Tuple[float, float]]
            x,y coordinates defined in npoints, points expected to be ordered clockwise by default

        filename : str, optional
            Name of avs polyline file to create, by default "polyline.inp"

        name : str, optional
            Internal lagrit name for mesh object, by default

        Returns
        -------
        MO

        Examples
        --------
        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()

        >>> mo = lg.tri_mo_from_polyline(
        ...     [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]
        ... )

        >>> lg.close()

        """
        mstr = str(len(coords)) + " " + str(len(coords)) + " 0 0 0\n"
        for i, p in enumerate(coords):
            mstr += " ".join([str(i + 1), str(p[0]), str(p[1]), str(0.0)])
            mstr += "\n"
        es1 = np.arange(len(coords)) + 1
        es2 = np.roll(es1, len(coords) - 1)
        for e1, e2 in zip(es1, es2):
            mstr += " ".join([str(e1), "1 line ", str(e1), str(e2)])
            mstr += "\n"
        with open(filename, "w") as fh:
            fh.write(mstr)
        # Check if name was specified, if not just generate one.
        if name is None:
            name = new_name("mo", self.core.mo_names())
        motmp = cast(MO, self.read_mo(filename))
        motri = motmp.copypts(elem_type="tri")
        motmp.delete()
        return motri

    def createpts(
        self,
        crd: Literal["xyz", "rtz", "rtp"],
        npts: Tuple[int, int, int],
        mins: Tuple[float, float, float],
        maxs: Tuple[float, float, float],
        elem_type: Literal[
            "tet", "hex", "pri", "pyr", "tri", "qua", "hyb", "lin", "triplane"
        ],
        vc_switch=(1, 1, 1),
        rz_switch=(1, 1, 1),
        rz_value=(1, 1, 1),
        connect=False,
        name: Optional[str] = None,
    ) -> "MO":
        """
        Create and Connect Points

        Parameters
        ----------
        crd : str
            Coordinate type of either 'xyz', 'rtz', or 'rtp'.
            Options include:
                xyz: Cartesian coordinates
                rtz: Cylindrical coordinates
                rtp: Spherical coordinates

        npts : Tuple[int, int, int]
            Number of points in each dimension.

        mins : Tuple[float, float, float]
            The starting values for each dimension.

        maxs : Tuple[float, float, float]
            The ending values for each dimension.

        elem_type : str
            The type of mesh object to create.

        vc_switch : Tuple[int, int, int], optional
            Determines if nodes represent vertices (1) or cell centers (0), by default (1, 1, 1)

        rz_switch : Tuple[int, int, int], optional
            Determines if ratio zoning is used (1) or not (0), by default (1, 1, 1)

        rz_value : Tuple[int, int, int], optional
            Ratio zoning values, by default (1, 1, 1)

        connect : bool, optional
            Whether or not to connect points, by default False

        name : str, optional
            Name of the mesh object, by default None

        Returns
        -------
        MO

        """
        if elem_type.startswith(("triplane", "qua")):
            assert (
                np.where(np.array(npts) <= 1)[0].shape[0] == 1
            ), f"{elem_type} elem_type requires one (1) in npts"  # noqa: S101
            assert (  # noqa: S101
                np.where((np.array(maxs) - np.array(mins)) == 0)[0][0] == 1
            ), f"{elem_type} elem_type requires one zero range (max-min)"
        if elem_type.startswith(("tet", "pri", "pyr", "hex")):
            assert np.all(
                np.array(npts) > 1
            ), f"{elem_type} elem_type requires all npts greater than 1"  # noqa: S101
            assert np.all((np.array(maxs) - np.array(mins)) > 0), (  # noqa: S101
                f"{elem_type} elem_type requires all ranges (max-min) greater than 0"
            )
        mo = self.create(elem_type=elem_type, name=name)
        mo.createpts(
            crd,
            npts,
            mins,
            maxs,
            vc_switch=vc_switch,
            rz_switch=rz_switch,
            rz_value=rz_value,
            connect=connect,
        )
        return mo

    def createpts_dxyz(
        self,
        dxyz: Tuple[float, float, float],
        mins: Tuple[float, float, float],
        maxs: Tuple[float, float, float],
        elem_type: Literal[
            "tet", "hex", "pri", "pyr", "tri", "qua", "hyb", "lin", "triplane"
        ],
        clip: Literal["under", "over"]
        | Tuple[
            Literal["under", "over"], Literal["under", "over"], Literal["under", "over"]
        ] = "under",
        hard_bound: Literal["min", "max"]
        | Tuple[
            Literal["min", "max"], Literal["min", "max"], Literal["min", "max"]
        ] = "min",
        rz_switch=(1, 1, 1),
        rz_value=(1, 1, 1),
        connect=True,
        name: Optional[str] = None,
    ):
        """
        Create and Connect Points to create an orthogonal hexahedral mesh. The
        vertex spacing is based on dxyz and the mins and maxs specified. mins
        (default, see hard_bound option) or maxs will be adhered to, while maxs
        (default) or mins will be modified based on the clip option to be
        truncated at the nearest value 'under' (default) or 'over' the range
        maxs-mins. clip and hard_bound options can be mixed by specifying tuples
        (see description below).

        Parameters
        ----------
        dxyz : Tuple[float, float, float]
            The spacing between points in x, y, and z directions.

        mins : Tuple[float, float, float]
            The starting value for each dimension.

        maxs : Tuple[float, float, float]
            The ending value for each dimension.

        elem_type : str
            The type of mesh object to create, automatically set to 'triplane' if 2d or 'tet' if 3d.

        clip : str or Tuple[str, str, str], optional
            How to handle bounds if range does not divide by dxyz, either 'under' or 'over' range, by default 'under'

        hard_bound : str or Tuple[str, str, str], optional

        rz_switch : Tuple[int, int, int], optional
            Determines if ratio zoning is used (1) or not (0), by default (1, 1, 1)

        rzvalue : Tuple[int, int, int], optional
            Ratio zoning values, by default (1, 1, 1)

        connect : bool, optional
            Whether or not to connect points, by default True

        name : str, optional
            Name of the mesh object, by default None

        Examples
        --------
        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()

        >>> # Create 2x2x2 cell mesh
        >>> mo = lg.create()
        >>> mo.createpts_dxyz(
        ...     (0.5, 0.5, 0.5),
        ...     (0.0, 0.0, 0.0),
        ...     (1.0, 1.0, 1.0),
        ...     rz_switch=[1, 1, 1],
        ...     connect=True,
        ... )

        >>> # Create 2x2x2 mesh where maxs will be truncated to nearest value under given maxs
        >>> mo_under = lg.create()
        >>> mo_under.createpts_dxyz(
        ...     (0.4, 0.4, 0.4),
        ...     (0.0, 0.0, 0.0),
        ...     (1.0, 1.0, 1.0),
        ...     rz_switch=[1, 1, 1],
        ...     connect=True,
        ... )

        >>> # Create 3x3x3 mesh where maxs will be truncated to nearest value over given maxs
        >>> mo_over = lg.create()
        >>> mo_over.createpts_dxyz(
        ...     (0.4, 0.4, 0.4),
        ...     (0.0, 0.0, 0.0),
        ...     (1.0, 1.0, 1.0),
        ...     clip="over",
        ...     rz_switch=[1, 1, 1],
        ...     connect=True,
        ... )

        >>> # Create 3x3x3 mesh where x and y maxs will be truncated to nearest value over given maxs
        >>> # and z min will be truncated  to nearest value
        >>> mo_mixed = lg.create()
        >>> mo_mixed.createpts_dxyz(
        ...     (0.4, 0.4, 0.4),
        ...     (0.0, 0.0, -1.0),
        ...     (1.0, 1.0, 0.0),
        ...     hard_bound=("min", "min", "max"),
        ...     clip=("under", "under", "over"),
        ...     rz_switch=[1, 1, 1],
        ...     connect=True,
        ... )

        >>> lg.close()

        """
        mo = self.create(elem_type=elem_type, name=name)
        mo.createpts_dxyz(
            dxyz,
            mins,
            maxs,
            clip=clip,
            hard_bound=hard_bound,
            rz_switch=rz_switch,
            rz_value=rz_value,
            connect=connect,
        )
        return mo

    def createpts_line(
        self,
        npts: Tuple[int, int, int],
        mins: Tuple[float, float, float],
        maxs: Tuple[float, float, float],
        vc_switch=(1, 1, 1),
        rz_switch=(1, 1, 1),
        name: Optional[str] = None,
    ) -> "MO":
        """
        Create and Connect Points in a line

        Parameters
        ----------

        npts : Tuple[int, int, int]
            The number of points in each dimension.

        mins : Tuple[float, float, float]
            The starting value for each dimension.

        maxs : Tuple[float, float, float]
            The ending value for each dimension.

        vc_switch : Tuple[int, int, int], optional
            Determines if nodes represent vertices (1) or cell centers (0), by default (1, 1, 1)

        rz_switch : Tuple[int, int, int], optional
            Determines if ratio zoning is used (1) or not (0), by default (1, 1, 1)

        Returns
        -------
        MO

        """
        mo = self.create("lin", name=name)
        mo.createpts_line(npts, mins, maxs, vc_switch=vc_switch, rz_switch=rz_switch)  # type: ignore

        return mo

    def gridder(
        self,
        x: Optional[List[float] | npt.NDArray[np.float64]] = None,
        y: Optional[List[float] | npt.NDArray[np.float64]] = None,
        z: Optional[List[float] | npt.NDArray[np.float64]] = None,
        connect=False,
        elem_type: Literal[
            "tet", "hex", "pri", "pyr", "tri", "qua", "hyb", "lin", "triplane"
        ] = "tet",
        name: Optional[str] = None,
        filename="gridder.inp",
    ):
        """
        Generate a logically rectangular orthogonal mesh corresponding to vectors of nodal positions.

        Parameters
        ----------
        x : List[float], optional
            x discretization locations

        y : List[float], optional
            y discretization locations

        z : List[float], optional
            z discretization locations

        connect : bool, optional
            Should the points be connected, by default False

        elem_type : str, optional
            Type of element for created mesh object, by default "tet"

        filename : str, optional
            Name of avs file created with nodal coordinates, by default "gridder.inp"

        Returns
        -------
        MO

        Examples
        --------
        >>> import pylagrit
        >>> import numpy as np
        >>> lg = pylagrit.PyLaGriT()

        >>> x0 = -np.logspace(1, 2, 15, endpoint=True)
        >>> x1 = np.arange(-10, 10, 1)
        >>> x2 = -x0
        >>> x = np.concatenate([x0, x1, x2])
        >>> y = x
        >>> mqua = lg.gridder(x, y, elem_type="qua", connect=True)

        >>> lg.close()

        """
        # TODO: validation for point set
        # dim = 0
        # if x is not None:
        #     if len(x) > 0:
        #         dim += 1
        # if y is not None:
        #     if len(y) > 0:
        #         dim += 1
        # if z is not None:
        #     if len(z) > 0:
        #         dim += 1
        # if dim == 0:
        #     print("ERROR: must define at least one of x, y, z arrays")
        #     return
        # if elem_type in ["line"] and dim != 1:
        #     print(
        #         "Error: Only 1 coordinate array (x,y,z) required for elem_type 'line'"
        #     )
        #     return
        # if elem_type in ["tri", "qua"] and dim != 2:
        #     print(
        #         "Error: Only 2 coordinate arrays (x,y,z) required for elem_type '"
        #         + str(elem_type)
        #         + "'"
        #     )
        #     return
        # if elem_type in ["tet", "hex"] and dim != 3:
        #     print(
        #         "Error: 3 coordinate arrays (x,y,z) required for elem_type '"
        #         + str(elem_type)
        #         + "'"
        #     )
        #     print("Set elem_type to a 2D format like 'qua' or 'triplane'")
        #     return
        xx = np.array([0]) if x is None or len(x) == 0 else np.unique(x)
        yy = np.array([0]) if y is None or len(y) == 0 else np.unique(y)
        zz = np.array([0]) if z is None or len(z) == 0 else np.unique(z)

        nodelist = np.array(list(product(*[zz, yy, xx])))
        nodelist = np.fliplr(nodelist)

        outfile = open(filename, "w")
        outfile.write(f"   {len(nodelist)} 0 0 0 0\n")
        for i, nd in enumerate(nodelist):
            outfile.write(f"{i:11d}        ")
            outfile.write(f"{nd[0]:14.8f}        ")
            outfile.write(f"{nd[1]:14.8f}        ")
            outfile.write(f"{nd[2]:14.8f}")
            outfile.write("\n")
        outfile.write("\n")
        outfile.close()

        m = (
            self.create(elem_type)
            if name is None
            else self.create(elem_type, name=name)
        )
        m.read(filename)

        if elem_type in ["qua", "hex"] and connect:
            cmd = [
                "createpts",
                "brick",
                "xyz",
                " ".join([str(len(xx)), str(len(yy)), str(len(zz))]),
                "1 0 0",
                "connect",
            ]
            m.sendcmd("/".join(cmd))
        elif connect:
            m.connect()

        self.sendcmd(f"cmo/printatt/{m.name}/-xyz- minmax")
        return m

    def points(
        self,
        coords: Points,
        connect=False,
        elem_type: Literal[
            "tet", "hex", "pri", "pyr", "tri", "qua", "hyb", "lin", "triplane"
        ] = "tet",
        filename="points.inp",
    ):
        """
        Generate a mesh object of points defined by x, y, z vectors.

        Parameters
        ----------
        coords : List[Tuple[float, float, float]]
            List of 3-tuples containing (x,y,z) coorinates

        connect : bool, optional
            Should the points be connected, by default False

        elem_type : str, optional
            Type of element for created mesh object, by default "tet"

        filename : str, optional
            Name of avs file created with nodal coordinates, by default "points.inp"

        Returns
        -------
        MO

        Examples
        --------
        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()

        >>> coords = [
        ...     [0, 0, 0],
        ...     [1, 0, 0],
        ...     [1, 1, 0],
        ...     [0, 1, 1],
        ...     [0, 0, 1],
        ...     [0, 1, 0],
        ...     [1, 1, 1],
        ...     [1, 0, 1],
        ... ]
        >>> m = lg.points(coords, elem_type="tet", connect=True)

        >>> lg.close()

        """
        # TODO: validation for point set
        # dim = 0
        # ix = numpy.all(numpy.diff(coords[:, 0]) == 0)
        # if not ix:
        #     dim += 1
        # iy = numpy.all(numpy.diff(coords[:, 1]) == 0)
        # if not iy:
        #     dim += 1
        # iz = numpy.all(numpy.diff(coords[:, 2]) == 0)
        # if not iz:
        #     dim += 1
        # if elem_type in ["line"] and dim != 1:
        #     print("Error: Coordinates must form line for elem_type 'line'")
        #     return
        # if elem_type in ["tri", "qua"] and dim != 2:
        #     print(
        #         "Error: Coordinates must form plane for elem_type '"
        #         + str(elem_type)
        #         + "'"
        #     )
        #     return
        # if elem_type in ["tet", "hex"] and dim != 3:
        #     print(
        #         "Error: 3D coordinates required for elem_type '" + str(elem_type) + "'"
        #     )
        #     print("Set elem_type to a 2D format like 'qua' or 'triplane'")
        #     return

        outfile = open(filename, "w")
        outfile.write("   " + str(len(coords)) + " 0 0 0 0\n")
        for i, nd in enumerate(coords):
            outfile.write(f"{i:11d}        ")
            outfile.write(f"{nd[0]:14.8f}        ")
            outfile.write(f"{nd[1]:14.8f}        ")
            outfile.write(f"{nd[2]:14.8f}")
            outfile.write("\n")
        outfile.write("\n")
        outfile.close()
        m = self.create(elem_type)
        m.read(filename)
        if elem_type in ["qua", "hex"] and connect:
            cmd = [
                "createpts",
                "brick",
                "xyz",
                " ".join([str(len(coords)), str(len(coords)), str(len(coords))]),
                "1 0 0",
                "connect",
            ]
            m.sendcmd("/".join(cmd))
        elif connect:
            m.connect()
        return m


class MO:
    """
    Mesh Object
    """

    def __init__(self, name: str, lg: PyLaGriT):
        self.lg = lg
        self.name = name
        self.obj = MeshObject(name)

        self.regions: Dict[str, Region] = {}
        self.mregions: Dict[str, MRegion] = {}
        self.surfaces: Dict[str, Surface] = {}

        lg.sendcmd("cmo/select/" + name)

    def select(self):
        self.lg.sendcmd("cmo/select/" + self.name)

    def __repr__(self):
        return self.name

    def sendcmd(self, cmd: str):
        self.select()
        self.lg.sendcmd(cmd)

    def attr(
        self, attr: str
    ) -> int | List[int] | float | List[float] | str | List[str]:
        return self.obj.attr(attr)

    @property
    def psets(self) -> List["PSet"]:
        return [PSet(k, self) for k in self.obj.pset_names()]

    @property
    def eltsets(self) -> List["EltSet"]:
        return [EltSet(k, self) for k in self.obj.eltset_names()]

    @property
    def elem_type(self):
        return self.obj.mesh_type()

    @property
    def points(self):
        x = self.attr("xic")
        y = self.attr("yic")
        z = self.attr("zic")

        return np.c_[x, y, z]

    def status(self, brief=False):
        cmd = ["cmo", "status", self.name]
        if brief:
            cmd.append("brief")

        self.sendcmd("/".join(cmd))

    def read(self, filename: str, filetype: Optional[str] = None):
        cmd = ["read"]

        if filetype is not None:
            cmd.append(filetype)

        # If filetype is lagrit, name is irrelevant
        if filetype != "lagrit":
            cmd.extend([filename, self.name])

        self.sendcmd("/".join(cmd))

    def printatt(
        self,
        attname: Optional[str] = None,
        stride=(1, 0, 0),
        pset: Optional["PSet"] = None,
        eltset: Optional["EltSet"] = None,
        ptype="value",
    ):
        stride = [str(v) for v in stride]

        if attname is None:
            attname = "-all-"

        if pset is not None:
            cmd = "/".join(
                [
                    "cmo/printatt",
                    self.name,
                    attname,
                    ptype,
                    ",".join(["pset", "get", str(pset)]),
                ]
            )
        elif eltset is not None:
            cmd = "/".join(
                [
                    "cmo/printatt",
                    self.name,
                    attname,
                    ptype,
                    ",".join(["eltset", "get", str(eltset)]),
                ]
            )
        else:
            cmd = "/".join(
                ["cmo/printatt", self.name, attname, ptype, ",".join(stride)]
            )

        self.sendcmd(cmd)

    def modatt(self, att: str, field: str, value: int | float | str):
        """
        Modify a field for a mesh object attribute
        """
        cmd = ["cmo", "modatt", str(self), att, field, str(value)]

        self.sendcmd("/".join(cmd))

    def delatt(self, attnames: Iterable[str], force=True):
        """
        Delete a list of attributes

        Parameters
        ----------
        attnames : Iterable[str]
            Attribute names to delete

        force : bool, optional
            If true, delete even if the attribute permanent persistance

        """
        # If single attribute as string, make list
        if isinstance(attnames, str):
            attnames = [attnames]
        for att in attnames:
            if force:
                cmd = "/".join(["cmo/DELATT", self.name, att])
            else:
                cmd = "/".join(["cmo/delatt", self.name, att])
            self.sendcmd(cmd)

    def copyatt(
        self,
        att_src: str,
        att_sink: Optional[str] = None,
        mo_src: Optional["MO"] = None,
    ):
        """
        Add a list of attributes

        Parameters
        ----------
        att_src: str
            Name of attribute to copy

        att_sink: str, optional
            Name of sink attribute

        mo_src: MO, optional
            Source mesh object

        """
        if att_sink is None:
            att_sink = att_src
        if mo_src is None:
            mo_src = self
        cmd = "/".join(["cmo", "copyatt", self.name, mo_src.name, att_sink, att_src])
        self.sendcmd(cmd)

    def addatt(
        self,
        attname: str,
        keyword: Optional[str] = None,
        vtype: Literal["vdouble", "vint"] = "vdouble",
        rank="scalar",
        length: Literal["nnodes", "nelements"] = "nnodes",
        interpolate="linear",
        persistence="permanent",
        ioflag="",
        value=0.0,
    ):
        """
        Add a list of attributes

        Parameters
        ----------
        attname: str
            Name of attribute to add

        keyword: str, optional
            Keyword used by lagrit for specific attributes

        vtype: str, optional
            Type of variable

        rank: str, optional

        length: str, optional
            Length of attribute

        """
        if keyword is not None:
            cmd = "/".join(["cmo/addatt", self.name, keyword, attname])
        else:
            cmd = "/".join(
                [
                    "cmo/addatt",
                    self.name,
                    attname,
                    vtype,
                    rank,
                    length,
                    interpolate,
                    persistence,
                    ioflag,
                    str(value),
                ]
            )
        self.sendcmd(cmd)

    def addatt_voronoi_volume(self, name="voronoi_volume"):
        """
        Add voronoi volume attribute to mesh object

        Parameters
        ----------
        name : str, optional
            Name of attribute to add, by default "voronoi_volume"

        """
        self.addatt(name, keyword="voronoi_volume")

    def addatt_voronoi_varea(self, attr_names="xvarea yvarea zvarea"):
        """
        Add voronoi area x,y,z component attributes for 2D planar mesh

        Parameters
        ----------
        attr_names : str
            Name of x,y,z attributes in LaGriT

        """
        self.addatt(attr_names, keyword="voronoi_varea")

    def minmax(self, attname: Optional[str] = None, stride=(1, 0, 0)):
        self.printatt(attname=attname, stride=stride, ptype="minmax")

    def minmax_xyz(self):
        cmd = "/".join(["cmo/printatt", self.name, "-xyz-", "minmax"])
        self.sendcmd(cmd)

    def list(
        self,
        attname: Optional[str] = None,
        stride=(1, 0, 0),
        pset: Optional["PSet"] = None,
    ):
        self.printatt(attname=attname, stride=stride, pset=pset, ptype="list")

    def setatt(self, attname: str, value: int | float, stride=(1, 0, 0)):
        stride = [str(v) for v in stride]
        cmd = "/".join(["cmo/setatt", self.name, attname, ",".join(stride), str(value)])
        self.sendcmd(cmd)

    def set_id(
        self,
        option: Literal["node", "element", "both"],
        node_attname="id_node",
        elem_attname="id_elem",
    ):
        """
        This command creates integer attributes that contain the node and/or
        element number. If later operations delete nodes or
        elements causing renumbering, these attributes will contain the
        original node or element number.

        Parameters
        ----------
        option : str
            create attribute for nodes, elements, or both

        node_attname : str
            name for new node attribute

        elem_attname : str
            name for new element attribute

        Examples
        -------
        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()

        Create source mesh
        >>> npts = (11, 11, 11)
        >>> mins = (0., 0., 0.)
        >>> maxs = (1., 1., 1.)
        >>> mesh = lg.create()
        >>> mesh.createpts_brick("xyz", npts, mins, maxs)

        Write node and element attribute numbers
        >>> mesh.set_id('both',node_attname='node_att1',elem_attname='elem_att1')

        >>> # select and remove points
        >>> p_mins = (0.5,0.,0.)
        >>> p_maxs = (1.,1.,1.)
        >>> points = mesh.pset_geom("xyz", p_mins, p_maxs)
        >>> mesh.rmpoint_pset(points)

        >>> lg.close()

        """
        if option == "both":
            cmd = "/".join(
                ["cmo/set_id", self.name, option, node_attname, elem_attname]
            )
        elif option == "node":
            cmd = "/".join(["cmo/set_id", self.name, option, node_attname])
        elif option == "element":
            cmd = "/".join(["cmo/set_id", self.name, option, elem_attname])

        self.sendcmd(cmd)

    def information(self):
        """
        Returns a formatted dictionary with mesh information.

        Information is that found in cmo/status/MO
        """

        atts = {
            "nodes": self.attr("nnodes"),
            "elmeents": self.attr("nelements"),
            "dimensions": self.attr("ndimensions_geom"),
            "type": self.elem_type[1],
            "dimensions_topology": self.attr("ndimensions_topo"),
            "attributes": {},
        }

        attr_data = self.obj.attr_list()

        for attr in attr_data:
            name = attr[0]
            atts["attributes"][name] = {}
            atts["attributes"][name]["type"] = attr[1]
            atts["attributes"][name]["rank"] = attr[2]
            atts["attributes"][name]["length"] = attr[3]
            atts["attributes"][name]["inter"] = attr[4]
            atts["attributes"][name]["persi"] = attr[5]
            atts["attributes"][name]["io"] = attr[6]
            # TODO: implement default value field for attribute
            atts["attributes"][name]["value"] = None

        return atts

    def pset_geom(
        self,
        mins: Tuple[float, float, float],
        maxs: Tuple[float, float, float],
        ctr: Tuple[float, float, float] = (0, 0, 0),
        geom: Literal["xyz", "rtz", "rtp"] = "xyz",
        stride: Tuple[int, int, int] = (1, 0, 0),
        name: Optional[str] = None,
    ):
        """
        Define PSet by Geometry

        Selects points from geometry specified by string geom and returns a
        PSet.

        Parameters
        ----------
        mins : Tuple[float, float, float]
            Coordinate of one of the shape's defining points.
                xyz (Cartesian):   (x1, y1, z1);
                rtz (Cylindrical): (radius1, theta1, z1);
                rtp (Spherical):   (radius1, theta1, phi1);

        maxs : Tuple[float, float, float]
            Coordinate of one of the shape's defining points.
                xyz (Cartesian):   (x2, y2, z2);
                rtz (Cylindrical): (radius2, theta2, z2);
                rtp (Spherical):   (radius2, theta2, phi2);

        ctr : Tuple[float, float, float], optional
            Coordinate of the relative center, by default (0, 0, 0)

        geom : str, optional
            Type of geometric shape, by default 'xyz'
                xyz: (spherical)
                rtz: (cylindrical)
                rtp: (spherical)

        stride : Tuple[int, int, int], optional
            Nodes defined by ifirst, ilast, and istride, by default (1, 0, 0)

        name : str, optional
            Name of the PSet, by default None

        Returns
        -------
        PSet

        """

        if name is None:
            name = new_name("p", self.lg.core.pset_names())

        cmd = "/".join(
            [
                "pset",
                name,
                "geom",
                geom,
                ",".join([str(v) for v in stride]),
                ",".join([str(v) for v in mins]),
                ",".join([str(v) for v in maxs]),
                ",".join([str(v) for v in ctr]),
            ]
        )
        self.sendcmd(cmd)
        return PSet(name, self)

    def pset_attribute(
        self,
        attribute: str,
        value: int | float,
        comparison: Literal["lt", "le", "gt", "ge", "eq", "ne"] = "eq",
        stride=(1, 0, 0),
        name: Optional[str] = None,
    ):
        """
        Define PSet by attribute

        Parameters
        ----------
        attribute : str
            Attribute name

        value : int | float
            Attribute value

        comparison : str, optional
            Attribute comparison, default is eq, by default "eq"

        stride : Tuple[int, int, int], optional
            Nodes defined by ifirst, ilast, and istride, by default (1, 0, 0)

        name : str, optional
            Name of the PSet created

        Returns
        -------
        PSet

        """
        if name is None:
            name = new_name("p", self.lg.core.pset_names())

        stride = [str(v) for v in stride]

        cmd = "/".join(
            [
                "pset",
                name,
                "attribute",
                attribute,
                ",".join(stride),
                str(value),
                comparison,
            ]
        )
        self.sendcmd(cmd)
        return PSet(name, self)

    def compute_distance(
        self,
        mo: "MO",
        option: Literal["distance_field", "signed_distance_field"] = "distance_field",
        attname="dfield",
    ):
        """
        Compute distance from one mesh object to another

        Parameters
        ----------
        mo : MO
            Mesh object to compute distance to base mesh from

        option : str, optional
            The type of distance field calculation

        attnames : str, optional
            Name of the attribute to be created in the base mesh, by default "dfield"

        Examples
        --------
        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()

        Create source mesh
        >>> npts = (1,91,1)
        >>> mins = (3.,0.,0.)
        >>> maxs = (3.,270.,0.)
        >>> src_mo = lg.create()
        >>> src_mo.createpts("rtz", npts, mins, maxs, connect=False)

        Create sink mesh
        >>> snk_mo = lg.create()
        >>> snk_mo.createpts("xyz", [30, 30, 1],[-5., -5., -5.],[5., 5., 5.], connect=False)

        Compute distance and store in sink mesh attribute 'dfield'
        >>> snk_mo.compute_distance(src_mo)
        >>> snk_mo.dump('comptest.gmv')

        >>> lg.close()

        """
        self.sendcmd("/".join(["compute", option, self.name, mo.name, attname]))

    def compute_extrapolate(
        self,
        surf_mo: "MO",
        dir: Literal["zpos", "zneg", "ypos", "yneg", "xpos", "xneg"] = "zpos",
        attname="zic",
    ):
        """
        Given a 3D mesh and a 2D surface, this command will extrapolate a scalar
        value from that surface onto every point of the mesh.

        Parameters
        ----------

        suf_mo : MO
            Surface mesh object to extrapolate from

        dir : str, optional
            The direction values are extrapolated from.

        attname : str, optional
            The name of the attribute in the surface mesh to be extrapolated

        Examples
        --------
        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()

        Create surface mesh
        >>> p1 = (-1.,-1.,-1.)
        >>> p2 = (301.,-1.,-1.)
        >>> p3 = (301.,301.,-1.)
        >>> p4 = (-1.,301.,-1.)
        >>> pts = [p1,p2,p3,p4]
        >>> nnodes = (30,30,1)
        >>> surf = lg.create("qua")
        >>> surf.quadxy(nnodes,pts)

        Make surface mesh interesting
        >>> surf.math('sin','zic',cmosrc=surf,attsrc='xic')
        >>> surf.math('multiply','zic',value=5.0,cmosrc=surf,attsrc='zic')
        >>> surf.perturb(0.,0.,1.)
        >>> surf.math('add','zic',value=60.0,cmosrc=surf,attsrc='zic')

        Create base mesh
        >>> hex = lg.create("hex")
        >>> hex.createpts_brick("xyz", (30, 30, 20), (0., 0., 0.), (300., 300., 50.))
        >>> hex.resetpts_itp()

        Extrapolate z values from surface mesh to base mesh
        >>> hex.compute_extrapolate(surf)
        >>> hex.dump('extrapolated.gmv')

        >>> lg.close()

        """
        self.sendcmd(
            "/".join(
                ["compute", "linear_transform", self.name, surf_mo.name, dir, attname]
            )
        )

    def pset_region(
        self, region: "Region", stride=(1, 0, 0), name: Optional[str] = None
    ):
        """
        Define PSet by region

        Parameters
        ----------
        region : Region
            Region to create pset

        stride : Tuple[int, int, int], optional
            Nodes defined by ifirst, ilast, and istride, by default (1, 0, 0)

        name : str, optional
            Name of the PSet created

        Returns
        -------
        PSet

        """
        if name is None:
            name = new_name("p", self.lg.core.pset_names())

        stride = [str(v) for v in stride]

        cmd = "/".join(["pset", name, "region", region.name, ",".join(stride)])
        self.sendcmd(cmd)
        return PSet(name, self)

    def pset_surface(
        self, surface: "Surface", stride=(1, 0, 0), name: Optional[str] = None
    ):
        """
        Define PSet by surface

        Parameters
        ----------
        surface : Surface
            Surface to create pset

        stride : Tuple[int, int, int], optional
            Nodes defined by ifirst, ilast, and istride, by default (1, 0, 0)

        name : str, optional
            Name of the PSet created

        Returns
        -------
        PSet

        """
        if name is None:
            name = new_name("p", self.lg.core.pset_names())

        stride = [str(v) for v in stride]

        cmd = "/".join(["pset", name, "surface", surface.name, ",".join(stride)])
        self.sendcmd(cmd)
        return PSet(name, self)

    def pset_bool(
        self,
        pset_list: Iterable["PSet"],
        boolean: Literal["union", "inter", "not"] = "union",
        name: Optional[str] = None,
    ):
        """
        Return PSet from boolean operation on list of psets

        Defines and returns a PSet from points that are not inside the PSet, ps.
        """
        # Generated a name if one is not specified.
        if name is None:
            name = new_name("p", self.lg.core.pset_names())

        # Create the new PSET in lagrit and the pylagrit object.
        cmd = ["pset", name, boolean]
        cmd.append(",".join([p.name for p in pset_list]))
        self.sendcmd("/".join(cmd))
        return PSet(name, self)

    def resetpts_itp(self):
        """
        set node type from connectivity of mesh

        """
        self.sendcmd("resetpts/itp")

    def eltset_object(
        self, mo: "MO", name: Optional[str] = None, attr_name: Optional[str] = None
    ):
        """
        Create element set from the intersecting elements with another mesh object
        """
        if name is None:
            name = new_name("e", self.lg.core.eltset_names())
        attr_name = self.intersect_elements(mo, attr_name)
        _ = self.eltset_attribute(attr_name, 0, boolstr="gt")
        return EltSet(name, self)

    def eltset_bool(
        self,
        eset_list: Iterable["EltSet"],
        boolstr: Literal["union", "inter", "not"] = "union",
        name: Optional[str] = None,
    ):
        """
        Create element set from boolean operation of set of element sets

        Parameters
        ----------
        eset_list : Iterable[EltSet]
            List of element sets to perform boolean operation

        boolstr : str
            Type of boolean operation to perform on element sets

        name : str, optional
            Name of the element set

        Returns
        -------
        EltSet

        """
        if name is None:
            name = new_name("e", self.lg.core.eltset_names())
        cmd = ["eltset", name, boolstr, " ".join([e.name for e in eset_list])]
        self.sendcmd("/".join(cmd))
        return EltSet(name, self)

    def eltset_region(self, region: "Region", name: Optional[str] = None):
        if name is None:
            name = new_name("e", self.lg.core.eltset_names())
        cmd = "/".join(["eltset", name, "region", region.name])
        self.sendcmd(cmd)
        return EltSet(name, self)

    def eltset_attribute(
        self,
        attribute_name: str,
        attribute_value: int | float,
        boolstr="eq",
        name: Optional[str] = None,
    ):
        if name is None:
            name = new_name("e", self.lg.core.eltset_names())
        cmd = "/".join(["eltset", name, attribute_name, boolstr, str(attribute_value)])
        self.sendcmd(cmd)
        return EltSet(name, self)

    def eltset_write(
        self, filename_root: str, eltset: Optional["EltSet"] = None, ascii=True
    ):
        """
        Write element set(s) to a file in ascii or binary format

        Parameters
        ----------
        filename_root : str
            Root name of file to write

        eltset : EltSet, optional
            Eltset to write; if blank, all eltsets in mesh object are written

        ascii : bool, optional
            Switch to indicate ascii (True) or binary (False) format, by default True

        Examples
        --------
        >>> import pylagrit
        >>> import numpy as np
        >>> lg = pylagrit.PyLaGriT()

        >>> dxyz = np.array([0.1, 0.25, 0.25])
        >>> mins = np.array([0.0, 0.0, 0.0])
        >>> maxs = np.array([1.0, 1.0, 1.0])
        >>> mqua = lg.createpts_dxyz(
        ...     dxyz, mins, maxs, "qua", hard_bound=("min", "max", "min"), connect=True
        ... )

        >>> example_pset1 = mqua.pset_geom("xyz", mins, maxs - (maxs - mins) / 2)
        >>> example_eset1 = example_pset1.eltset()
        >>> example_pset2 = mqua.pset_geom("xyz", mins + maxs / 2, maxs)
        >>> example_eset2 = example_pset2.eltset()

        Write specific eltset
        >>> mqua.eltset_write("test_specific", eltset=example_eset1)

        Write all eltsets
        >>> mqua.eltset_write("test_all")

        >>> lg.close()

        """
        name = "-all-" if eltset is None else eltset.name
        ascii = "ascii" if ascii else "binary"
        cmd = "/".join(["eltset", name, "write", filename_root, ascii])
        self.lg.sendcmd(cmd)

    def rmpoint_pset(
        self, pset: "PSet", itype="exclusive", compress=True, resetpts_itp=True
    ):
        cmd = "rmpoint/pset,get," + pset.name + "/" + itype
        self.sendcmd(cmd)
        if compress:
            self.rmpoint_compress(resetpts_itp=resetpts_itp)

    def rmpoint_eltset(self, eltset: "EltSet", compress=True, resetpts_itp=True):
        cmd = "rmpoint/element/eltset,get," + eltset.name
        self.sendcmd(cmd)
        if compress:
            self.rmpoint_compress(resetpts_itp=resetpts_itp)

    def rmpoint_compress(self, filter_bool=False, resetpts_itp=True):
        """
        Remove all marked nodes and correct the itet array

        Parameters
        ----------
        resetpts_itp : bool
            set node type from connectivity of mesh

        """

        if filter_bool:
            self.sendcmd("filter/1,0,0")
        self.sendcmd("rmpoint/compress")
        if resetpts_itp:
            self.resetpts_itp()

    def reorder_nodes(self, order="ascending", cycle="zic yic xic"):
        self.sendcmd("resetpts itp")
        self.sendcmd("/".join(["sort", self.name, "index", order, "ikey", cycle]))
        self.sendcmd("reorder / " + self.name + " / ikey")
        self.sendcmd("cmo / DELATT / " + self.name + " / ikey")

    def trans(
        self,
        xold: Tuple[float, float, float],
        xnew: Tuple[float, float, float],
        stride=(1, 0, 0),
    ):
        """Translate mesh according to old coordinates "xold" to new coordinates "xnew"

        Parameters
        ----------
        xold : Tuple[float, float, float]
            old position

        xnew : Tuple[float, float, float]
            new position

        stride : Tuple[int, int, int], optional
            Nodes defined by ifirst, ilast, and istride, by default (1, 0, 0)

        """
        cmd = "/".join(
            [
                "trans",
                ",".join([str(v) for v in stride]),
                ",".join([str(v) for v in xold]),
                ",".join([str(v) for v in xnew]),
            ]
        )
        self.sendcmd(cmd)

    def rotateln(
        self,
        coord1: Tuple[float, float, float],
        coord2: Tuple[float, float, float],
        theta: float,
        center=(0.0, 0.0, 0.0),
        copy=False,
        stride=(1, 0, 0),
    ):
        """
        Rotates a point distribution (specified by ifirst,ilast,istride) about a line.
        The copy option allows the user to make a copy of the original points as well
        as the rotated points, while copy=False just keeps the rotated points themselves.
        The line of rotation defined by coord1 and coord2 needs to be defined such that
        the endpoints extend beyond the point distribution being rotated. theta (in degrees)
        is the angle of rotation whose positive direction is determined by the right-hand-rule,
        that is, if the thumb of your right hand points in the direction of the line
        (1 to 2), then your fingers will curl in the direction of rotation. center is the point
        where the line can be shifted to before rotation takes place.
        If the copy option is chosen, the new points will have only coordinate values
        (xic, yic, zic); no values will be set for any other mesh object attribute for these points.
        Note:  The end points of the  line segment must extend well beyond the point set being rotated.

        Examples
        --------
        Example 1

        >>> import pylagrit
        >>> import numpy as np
        >>> lg = pylagrit.PyLaGriT()

        >>> x = np.arange(0, 10.1, 1)
        >>> y = x
        >>> z = [0, 1]

        >>> mqua = lg.gridder(x, y, z, elem_type="hex", connect=True)
        >>> mqua.rotateln((mqua.attr("xmin") - 0.1, 0, 0), (mqua.attr("xmax") + 0.1, 0, 0), 25)
        >>> mqua.dump_exo("rotated.exo")
        >>> mqua.dump_ats_xml("rotated.xml", "rotated.exo")

        >>> lg.close()

        Example 2

        >>> import pylagrit
        >>> import numpy as np
        >>> lg = pylagrit.PyLaGriT()

        >>> x = np.arange(0, 10.1, 1)
        >>> y = [0, 1]

        >>> layer = lg.gridder(x=x, y=y, elem_type="qua", connect=True)
        >>> layer.rotateln([0, layer.attr("ymin") - 0.10, 0], [0, layer.attr("ymax") + 0.1, 0], 25)
        >>> layer.dump("tmp_lay_top.inp")

        >>> layers = [0.1, 1.0]
        >>> addnum = [4, 2]
        >>> matnum = [2, 1]
        >>> layer_interfaces = np.cumsum(layers)
        >>> mtop = layer.copy()
        >>> stack_files = ["tmp_lay_top.inp 1,9"]

        >>> i = 1
        >>> for li, m, a in zip(layer_interfaces, matnum, addnum):
        ...     layer.math("sub", "zic", li, cmosrc=mtop)
        ...     stack_files.append(f"tmp_lay{str(i)}.inp {str(m)}, {str(a)}")
        ...     layer.dump(f"tmp_lay{str(i)}.inp")
        ...     i += 1
        >>> layer.math("sub", "zic", 2, cmosrc=mtop)

        >>> layer.dump("tmp_lay_bot.inp")
        >>> stack_files.append("tmp_lay_bot.inp 2")
        >>> stack_files.reverse()

        Create stacked layer mesh and fill
        >>> stack_hex = lg.create()
        >>> stack_hex.stack_layers(stack_files, "avs", flip_opt=True)
        >>> stack = stack_hex.stack_fill()

        # TODO: Exodus suppport
        # >>> stack.dump_exo("rotated.exo")
        # >>> stack.dump_ats_xml("rotated.xml", "rotated.exo")

        >>> lg.close()

        """
        self.sendcmd(
            "/".join(
                [
                    "rotateln",
                    ",".join([str(v) for v in stride]),
                    "copy" if copy else "nocopy",
                    ",".join([str(v) for v in coord1]),
                    ",".join([str(v) for v in coord2]),
                    str(theta),
                    ",".join([str(v) for v in center]),
                ]
            )
        )

    def massage(
        self,
        bisection_len: float,
        merge_len: float,
        toldamage: float,
        tolroughness: Optional[float] = None,
        stride: Optional[Tuple[int, int, int]] = None,
        nosmooth=False,
        norecon=False,
        strictmergelength=False,
        checkaxy=False,
        semiexclusive=False,
        ignoremats=False,
        lite=False,
    ):
        """
        MASSAGE creates, annihilates, and moves nodes and swaps connections in a 2D or 3D mesh
        in order to improve element aspect ratios and establish user-desired edge lengths.

        The actions of MASSAGE are controlled by values of these four parameters:

            bisection_length  - edge length that will trigger bisection.
            merge_length - edge length that will trigger merging.
            toldamage - maximum grid deformation of interfaces and external boundaries
                        allowed in a single merge, smooth or reconnection event.
            tolroughness - (for 2D surface grids only)  measure of grid roughness
                           (deviation from average surface normal) that triggers refinement.

        The final, optional keywork argument(s) can be one or more of nosmooth, norecon, lite,
        ignoremats, strictmergelength, checkaxy, semiexclusive, and exclusive.

        Specifying nosmooth will turn off the 'smooth' step by skipping the call to SGD.
        Specifying norecon will turn off all 'recon' steps.
        If lite is specified, only one iteration of the merging/reconnection/smoothing
        loop is executed, and a reconnection after edge refinement is omitted.
        This is suitable for applications, such as Gradient Weighted Moving Finite
        Elements, where MASSAGE is called repeatedly.

        The optional argument ignoremats causes MASSAGE to process the multimaterial
        mesh in a single material mode; it ignores the material interfaces.

        The optional argument strictmergelength forces strict interpretation of
        merge_length so that there is no merging along the edges of flat elements.
        This is important if ignoremats is specified to avoid losing the interfaces.

        If checkaxy is given, then we insure that for 2D meshes, the output mesh
        will have positive xy-projected triangle areas, provided that the input mesh
        had them in the first place.

        If exclusive is given, then edge refinement operations will only be performed
        on edges whose endpoints are both in the PSET that MASSAGE is working on.
        (As usual, new nodes created by refinement are added to the PSET so that MASSAGE
        can refine edges recursively.)  The default behavior is 'inclusive',
        where only ONE edge endpoint has to belong to the PSET for the edge to be
        eligible for refinement.

        If semiexclusive is given, refinement will only be triggered by edges with
        both endpoints in the PSET, but some edges with less than two endpoints in
        the PSET might be refined as part of a 'Rivara chain' triggered by the refinement
        of an edge with both endpoints in the PSET.  This represents an intermediate
        case between 'inclusive' and exclusive
        """

        cmd = ["massage", str(bisection_len), str(merge_len), str(toldamage)]

        if tolroughness is not None:
            cmd.append(str(tolroughness))
        if stride is not None:
            cmd.append(",".join([str(x) for x in stride]))

        # Add optional boolean arguments
        _iter = zip(
            [
                "nosmooth",
                "norecon",
                "strictmergelength",
                "checkaxy",
                "semiexclusive",
                "ignoremats",
                "lite",
            ],
            [
                nosmooth,
                norecon,
                strictmergelength,
                checkaxy,
                semiexclusive,
                ignoremats,
                lite,
            ],
        )
        [cmd.append(c[0]) for c in _iter if c[1]]
        self.sendcmd("/".join(cmd))

    def massage2(
        self,
        filename: str,
        min_scale: float,
        bisection_len: float,
        merge_len: float,
        toldamage: float,
        tolroughness: Optional[float] = None,
        stride: Optional[Tuple[int, int, int]] = None,
        nosmooth=False,
        norecon=False,
        strictmergelength=False,
        checkaxy=False,
        semiexclusive=False,
        ignoremats=False,
        lite=False,
    ):
        """
        MASSAGE2 iteratively calls MASSAGE to refine adaptively according to a
        gradient field. Thus, the bisection_length option must be a field.

        file_name is a file which contains a set of LaGriT commands that
        calculates the gradient field based on the distance field. In other
        words, the gradient field is a function of the distance field.
        It is necessary to have this file when using this routine, as the field
        must be updated after each refinement iteration.

        Use this function in conjunction with PyLaGriT.define(**kwargs) for
        best results.

        See MASSAGE for other arguments.
        """

        cmd = [
            "massage2",
            filename,
            str(min_scale),
            str(bisection_len),
            str(merge_len),
            str(toldamage),
        ]
        if tolroughness is not None:
            cmd.append(str(tolroughness))
        if stride is not None:
            cmd.append(",".join([str(x) for x in stride]))

        # Add optional boolean arguments
        _iter = zip(
            [
                "nosmooth",
                "norecon",
                "strictmergelength",
                "checkaxy",
                "semiexclusive",
                "ignoremats",
                "lite",
            ],
            [
                nosmooth,
                norecon,
                strictmergelength,
                checkaxy,
                semiexclusive,
                ignoremats,
                lite,
            ],
        )
        [cmd.append(c[0]) for c in _iter if c[1]]
        self.sendcmd("/".join(cmd))

    def perturb(self, xfactor: float, yfactor: float, zfactor: float, stride=(1, 0, 0)):
        """
        This command moves node coordinates in the following manner.

        Three pairs of random numbers between 0 and 1 are generated.
        These pairs refer to the x, y and z coordinates of the nodes respectively.
        The first random number of each pair is multiplied by the factor given in
        the command. The second random number is used to determine
        if the calculated offset is to be added or subtracted from the coordinate.
        """

        cmd = [
            "perturb",
            ",".join([str(x) for x in stride]),
            str(xfactor),
            str(yfactor),
            str(zfactor),
        ]
        self.sendcmd("/".join(cmd))

    def upscale(
        self,
        method: Literal["ariave", "geoave", "harave", "min", "max", "sum"],
        attsink: str,
        cmosrc: "MO",
        attsrc: Optional[str] = None,
        stride=(1, 0, 0),
        boundary_choice: Optional[Literal["single", "divide", "multiple"]] = None,
        keepatt=False,
        set_id=False,
    ):
        """
        The upscale command is used to interpolate attribute values from nodes of a fine source mesh to node
        attributes of a coarse sink mesh. The subroutine finds nodes of the fine source mesh within the Voronoi
        cell of every node in the coarser sink mesh. Nodes on cell boundaries are assigned to two or more sink
        nodes. Then the attributes of all the source nodes within a source node's cell are upscaled into a
        single value based on the chosen method. Mesh elements and connectivity are ignored and only node
        values are used to upscale values on to the sink mesh nodes.

        Parameters
        ----------
        method : str
            Type of upscaling: sum, min, max, and averages ariave, harave, geoave

        attsink : str
            attribute sink

        cmosrc : MO
            PyLaGriT mesh object source

        attsrc : str, optional
            attribute src, defaults to name of attsink

        stride : Tuple[int, int, int], optional
            tuple of (first, last, stride) of points, by default (1, 0, 0)

        boundary_choice : str, optional
            method of choice when source nodes are found on the boundary of multiple Voronoi volumes of sink nodes

        """
        if attsrc is None:
            attsrc = attsink
        cmd = [
            "upscale",
            method,
            self.name,
            attsink,
            ",".join([str(v) for v in stride]),
            cmosrc.name,
            attsrc,
        ]
        opts = []
        if boundary_choice is not None:
            opts.append(boundary_choice)
        if keepatt:
            opts.append("keepatt")
        if set_id:
            opts.append("set_id")
        if len(opts) > 0:
            cmd.append(" ".join(opts))
        self.sendcmd("/".join(cmd))

    @property
    def vtk_celltypes(self):
        el_type = np.array(self.attr("itettyp"))

        return VTK_CELL_TYPE_MAP[el_type - 1]

    @property
    def cells(self):
        el_offset = cast(List[int], self.attr("itetoff"))
        el_nodes = np.array(self.attr("itet")) - 1

        celltypes = self.vtk_celltypes
        n_nodes = np.array([CELL_NODE_NUMBERS[t] for t in celltypes])

        cells = np.empty(self.attr("nelements") + n_nodes.sum(), dtype=np.int64)

        el_idx = 0
        for offset, npts, typ in zip(el_offset, n_nodes, celltypes):
            cells[el_idx] = npts

            tmp = el_nodes[offset : offset + npts]
            if typ == CellType.WEDGE:  # prism type in lagrit
                tmp = tmp[[0, 2, 1, 3, 5, 4]]

            cells[el_idx + 1 : el_idx + 1 + npts] = tmp
            el_idx = el_idx + 1 + npts

        return cells

    def to_pyvista(self):
        """
        Convert mesh object to pyvista unstructured grid
        """
        mesh = pv.UnstructuredGrid(self.cells, self.vtk_celltypes, self.points)

        info = self.information()
        export_attrs = [
            attr
            for attr in info["attributes"].items()
            # write only avs attributes
            if "a" in attr[1]["io"] and attr[0] != "-def-"
        ]

        for attr in export_attrs:
            attr_name = cast(str, attr[0])
            length = attr[1]["length"]
            attr_data = cast(List[int | float], self.attr(attr_name))

            if length == "nnodes":
                mesh.point_data[attr_name] = attr_data
            elif length == "nelements":
                mesh.cell_data[attr_name] = attr_data

        return mesh

    def plot(self, show_axes=True, show_edges=True, show_bounds=True, **kwargs):
        self.to_pyvista().plot(
            show_axes=show_axes,
            show_edges=show_edges,
            show_bounds=show_bounds,
            **kwargs,
        )

    def dump(self, filename: Optional[str] = None, format: Optional[str] = None, *args):
        if filename is None and format is None:
            raise ValueError(
                "At least one of either filename or format option is required"
            )

        if filename is not None and format is not None:
            if format in ["fehm", "zone_outside", "zone_outside_minmax"]:
                filename = filename.split(".")[0]
            if format == "stor" and len(args) == 0:
                filename = filename.split(".")[0]
            cmd = "/".join(["dump", format, filename, self.name])
        elif format is not None:
            if format in ["avs", "avs2"]:
                filename = self.name + ".inp"
            elif format == "fehm":
                filename = self.name
            elif format == "gmv":
                filename = self.name + ".gmv"
            elif format == "tecplot":
                filename = self.name + ".plt"
            elif format == "lagrit":
                filename = self.name + ".lg"
            elif format == "exo":
                filename = self.name + ".exo"
            else:
                raise NotImplementedError("Unsupported format")
            cmd = "/".join(["dump", format, filename, self.name])
        elif filename is not None:
            cmd = "/".join(["dump", filename, self.name])
        else:
            cmd = "/".join(["dump", self.name + ".inp"])

        for arg in args:
            cmd = "/".join([cmd, str(arg)])

        self.sendcmd(cmd)

    def dump_ats_xml(
        self,
        filename: str,
        meshfilename: str,
        matnames: Dict[int, str] = {},
        facenames: Dict[int, str] = {},
    ):
        """
        Write ats style xml file with regions

        Parameters
        ----------
        filename : str
            Name of xml file to write

        meshfilename : str
            Name of exodus file to use in xml

        matnames : dict, optional
            Dictionary of region names keyed by exodus material number

        facenames : dict, optional
            Dictionary of faceset names keyed by exodus faceset number

        """
        main = ET.Element("ParameterList", {"name": "Main", "type": "ParameterList"})

        ET.SubElement(
            main,
            "Parameter",
            {"name": "Native Unstructured Input", "type": "bool", "value": "true"},
        )
        ET.SubElement(
            main,
            "Parameter",
            {"name": "grid_option", "type": "string", "value": "Unstructured"},
        )

        mesh = ET.SubElement(
            main, "ParameterList", {"name": "Mesh", "type": "ParameterList"}
        )
        ET.SubElement(
            mesh,
            "Parameter",
            {"isUsed": "true", "name": "Framework", "type": "string", "value": "MSTK"},
        )

        mesh1 = ET.SubElement(
            mesh, "ParameterList", {"name": "Read Mesh File", "type": "ParameterList"}
        )
        ET.SubElement(
            mesh1,
            "Parameter",
            {"name": "File", "type": "string", "value": meshfilename},
        )
        ET.SubElement(
            mesh1,
            "Parameter",
            {"name": "Format", "type": "string", "value": "Exodus II"},
        )

        mesh2 = ET.SubElement(
            mesh, "ParameterList", {"name": "Surface Mesh", "type": "ParameterList"}
        )
        ET.SubElement(
            mesh2,
            "Parameter",
            {"name": "surface sideset name", "type": "string", "value": "surface"},
        )
        mesh2a = ET.SubElement(
            mesh2, "ParameterList", {"name": "Expert", "type": "ParameterList"}
        )
        ET.SubElement(
            mesh2a,
            "Parameter",
            {"name": "Verify Mesh", "type": "bool", "value": "false"},
        )

        r = ET.SubElement(
            main, "ParameterList", {"name": "Regions", "type": "ParameterList"}
        )

        r1 = ET.SubElement(
            r,
            "ParameterList",
            {"name": "computational domain", "type": "ParameterList"},
        )
        l1 = ET.SubElement(
            r1, "ParameterList", {"name": "Region: Box", "type": "ParameterList"}
        )
        ET.SubElement(
            l1,
            "Parameter",
            {
                "name": "Low Coordinate",
                "type": "Array(double)",
                "value": "{-1.e20,-1.e20,-1.e20}",
            },
        )
        ET.SubElement(
            l1,
            "Parameter",
            {
                "name": "High Coordinate",
                "type": "Array(double)",
                "value": "{1.e20,1.e20,1.e20}",
            },
        )

        r2 = ET.SubElement(
            r, "ParameterList", {"name": "surface domain", "type": "ParameterList"}
        )
        l2 = ET.SubElement(
            r2, "ParameterList", {"name": "Region: Box", "type": "ParameterList"}
        )
        ET.SubElement(
            l2,
            "Parameter",
            {
                "name": "Low Coordinate",
                "type": "Array(double)",
                "value": "{-1.e20,-1.e20}",
            },
        )
        ET.SubElement(
            l2,
            "Parameter",
            {
                "name": "High Coordinate",
                "type": "Array(double)",
                "value": "{1.e20,1.e20}",
            },
        )

        rmat = []
        lmat = []
        for k, v in matnames.items():
            rmat.append(
                ET.SubElement(
                    r, "ParameterList", {"name": str(v), "type": "ParameterList"}
                )
            )
            lmat.append(
                ET.SubElement(
                    rmat[-1],
                    "ParameterList",
                    {"name": "Region: Labeled Set", "type": "ParameterList"},
                )
            )
            ET.SubElement(
                lmat[-1],
                "Parameter",
                {"name": "Label", "type": "string", "value": str(k)},
            )
            ET.SubElement(
                lmat[-1],
                "Parameter",
                {"name": "File", "type": "string", "value": meshfilename},
            )
            ET.SubElement(
                lmat[-1],
                "Parameter",
                {"name": "Format", "type": "string", "value": "Exodus II"},
            )
            ET.SubElement(
                lmat[-1],
                "Parameter",
                {"name": "Entity", "type": "string", "value": "Cell"},
            )

        rsurf = []
        lsurf = []
        for k, v in facenames.items():
            rsurf.append(
                ET.SubElement(
                    r, "ParameterList", {"name": str(v), "type": "ParameterList"}
                )
            )
            lsurf.append(
                ET.SubElement(
                    rsurf[-1],
                    "ParameterList",
                    {"name": "Region: Labeled Set", "type": "ParameterList"},
                )
            )
            ET.SubElement(
                lsurf[-1],
                "Parameter",
                {"name": "Label", "type": "string", "value": str(k)},
            )
            ET.SubElement(
                lsurf[-1],
                "Parameter",
                {"name": "File", "type": "string", "value": meshfilename},
            )
            ET.SubElement(
                lsurf[-1],
                "Parameter",
                {"name": "Format", "type": "string", "value": "Exodus II"},
            )
            ET.SubElement(
                lsurf[-1],
                "Parameter",
                {"name": "Entity", "type": "string", "value": "Face"},
            )

        m_str = ET.tostring(main)
        m_reparsed = minidom.parseString(m_str)  # noqa: S318
        with open(filename, "w") as f:
            f.write(m_reparsed.toprettyxml(indent="  "))

    def dump_exo(
        self,
        filename: str,
        psets=False,
        eltsets=False,
        facesets: Optional[Iterable["FaceSet"]] = None,
    ):
        """
        Dump exo file

        Parameters
        ----------
        filename : str
            Name of exo file

        psets : bool, optional
            Boolean indicating that exodus will only include psets, by default False

        eltsets : bool, optional
            Boolean indicating that exodus will only include element sets, by default False

        facesets : Iterable[FaceSet], optional
            Array of FaceSet objects

        Examples
        --------
        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()

        >>> mo = lg.create()
        >>> mo.createpts(
        ...     "xyz",
        ...     (3, 3, 3),
        ...     (0.0, 0.0, 0.0),
        ...     (1.0, 1.0, 1.0),
        ...     rz_switch=[1, 1, 1],
        ...     connect=True,
        ... )
        >>> fs = mo.create_boundary_facesets(base_name="faceset_bounds")

        # TODO: Support exodus file
        # >>> mo.dump_exo("cube.exo", facesets=fs.values())

        >>> lg.close()

        """
        cmd = ["dump", "exo", filename, self.name]
        if psets:
            cmd.append("psets")
        if eltsets:
            cmd.append("eltsets")
        if facesets is not None:
            cmd.append("facesets")
            cmd.extend([fc.filename for fc in facesets])

        self.sendcmd("/".join(cmd))

    def dump_zone_imt(self, filename: str, imt_value: int):
        cmd = ["dump", "zone_imt", filename, self.name, str(imt_value)]
        self.sendcmd("/".join(cmd))

    def dump_zone_outside(
        self, filename: str, keepatt=False, keepatt_median=False, keepatt_voronoi=False
    ):
        if keepatt_median and keepatt_voronoi:
            raise ValueError("keepatt_median and keepatt_voronoi cannot both be True")

        cmd = ["dump", "zone_outside", filename, self.name]
        if keepatt:
            cmd.append("keepatt")
        elif keepatt_median:
            cmd.append("keepatt_median")
        elif keepatt_voronoi:
            cmd.append("keepatt_voronoi")
        self.sendcmd("/".join(cmd))

    def dump_pset(
        self,
        filerootname: str,
        zonetype: Literal["zone", "zonn"] = "zone",
        pset: List["PSet"] = [],
    ):
        """
        Dump zone file of psets

        Parameters
        ----------
        filename_root : str
            root name of files to create, pset name will be added to name

        zonetype : str
            Type of zone file to dump, 'zone' or 'zonn'

        pset : List[PSet]
            list of psets to dump, all psets dumped if empty list

        """
        if len(pset) == 0:
            cmd = ["pset", "-all-", zonetype, filerootname, "ascii"]
            self.sendcmd("/".join(cmd))
        else:
            for p in pset:
                cmd = ["pset", p.name, zonetype, filerootname + "_" + p.name, "ascii"]
                self.sendcmd("/".join(cmd))

    def dump_pflotran_ugi(self, filename: str):
        """
        Dump PFLOTRAN unstructured implicit (ugi) file
        """
        points = self.points
        eltypeids = cast(int, self.attr("itettyp"))
        eltypes = PF_CELL_TYPE_MAP[eltypeids - 1]
        cells = self.cells

        with open(filename, "w") as pf_ugi:
            pf_ugi.write(f"{self.attr('nelements')} {self.attr('nnodes')}\n")

            cell_idx = 0
            for n in range(cast(int, self.attr("nelements"))):
                if eltypes[n] is None:
                    raise ValueError("unsupported element type:")

                pf_ugi.write(f"{eltypes[n]} ")

                n_nodes = cells[cell_idx]
                cell_idx += 1
                pf_ugi.write(" ".join(map(str, cells[cell_idx : cell_idx + n_nodes])))
                pf_ugi.write("\n")
                cell_idx += n_nodes

            for p in points:
                pf_ugi.write(" ".join(map(str, p)))
                pf_ugi.write("\n")

    def dump_pflotran_surface(self, surface_mesh: "MO", filename_prefix: str):
        """
        Dump PFLOTRAN implicit surface definition (*.ss) file and return the surface mesh object
        """
        el_ids = cast(List[int], surface_mesh.obj.attr("idelem1"))
        elface_ids = cast(List[int], surface_mesh.obj.attr("idface1"))
        face_ids = cast(List[int], surface_mesh.obj.attr("itetclr"))

        cells = self.cells
        stacked_cells = []
        cell_idx = 0
        for _ in range(cast(int, self.attr("nelements"))):
            n_nodes = cells[cell_idx]
            stacked_cells.append(
                # cell id starts from 1 in PFLOTRAN
                np.array(cells[cell_idx : cell_idx + 1 + n_nodes]) + 1
            )
            cell_idx += 1 + n_nodes

        eltypeids = cast(List[int], self.attr("itettyp"))

        face_collection = defaultdict(list)
        for el, elface, face in zip(el_ids, elface_ids, face_ids):
            eltype = eltypeids[el - 1]
            face_map = NODES_FOR_FACE[eltype - 1]
            if face_map is None:
                continue

            mask = face_map[elface - 1]
            face_collection[face].append(stacked_cells[el - 1][mask])

        for face_id, nodes in face_collection.items():
            with open(f"{filename_prefix}_{face_id:02}.ss", "w") as f:
                f.write(f"{len(nodes)}\n")
                for node_ids in nodes:
                    if len(node_ids) == 3:
                        f.write("T ")
                    elif len(node_ids) == 4:
                        f.write("Q ")
                    else:
                        raise ValueError("unsupported boundary faces")

                    f.write(" ".join(map(str, node_ids)))
                    f.write("\n")

    def dump_pflotran_uge(self, filename_root: str, nofilter_zero=False):
        """
        Dump PFLOTRAN UGE file

        Parameters
        ----------
        filename_root : str
            root name of UGE file

        nofilter_zero : bool
            Set to true to write zero coefficients to file

        Examples
        --------
        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()

        >>> mo = lg.create()
        >>> mo.createpts(
        ...     "xyz",
        ...     (3, 3, 3),
        ...     (0.0, 0.0, 0.0),
        ...     (1.0, 1.0, 1.0),
        ...     rz_switch=[1, 1, 1],
        ...     connect=True,
        ... )
        >>> mo.dump_pflotran_uge("test_pflotran_dump")

        >>> lg.close()

        """
        cmd = ["dump", "pflotran", filename_root, self.name]
        if nofilter_zero:
            cmd.append("nofilter_zero")
        self.sendcmd("/".join(cmd))

    def delete(self):
        self.sendcmd("cmo/delete/" + self.name)

    def create_boundary_facesets(
        self,
        stacked_layers=False,
        base_name: Optional[str] = None,
        reorder=False,
        external=True,
    ) -> OrderedDict[str, "FaceSet"]:
        """
        Creates facesets for each boundary and writes associated avs faceset file

        Parameters
        ----------
        base_name : str, optional
            base name of faceset files

        stacked_layers : bool, optional
            if mesh is created by stack_layers, user layertyp attr to determine top and bottom

        reorder : bool, optional
            reorder nodes on cell medians, usually needed for exodus file

        Returns
        -------
        OrderedDict[str, FaceSet]

        """
        if base_name is None:
            base_name = "faceset_" + self.name
        mo_surf = self.extract_surfmesh(reorder=reorder, external=external)
        mo_surf.addatt("id_side", vtype="vint", rank="scalar", length="nelements")
        mo_surf.settets("normal")
        mo_surf.copyatt("itetclr", "id_side")
        mo_surf.delatt(["id_side"])
        fs = OrderedDict()
        if stacked_layers:
            pbot = mo_surf.pset_attribute("layertyp", -1)
            ebot = pbot.eltset(membership="exclusive")
        else:
            ebot = mo_surf.eltset_attribute("itetclr", 1)
        fs["bottom"] = ebot.create_faceset(f"{base_name}_bottom.avs")
        if stacked_layers:
            ptop = mo_surf.pset_attribute("layertyp", -2)
            etop = ptop.eltset(membership="exclusive")
        else:
            etop = mo_surf.eltset_attribute("itetclr", 2)
        fs["top"] = etop.create_faceset(f"{base_name}_top.avs")
        er = mo_surf.eltset_attribute("itetclr", 3)
        fs["right"] = er.create_faceset(f"{base_name}_right.avs")
        eback = mo_surf.eltset_attribute("itetclr", 4)
        fs["back"] = eback.create_faceset(f"{base_name}_back.avs")
        el = mo_surf.eltset_attribute("itetclr", 5)
        fs["left"] = el.create_faceset(f"{base_name}_left.avs")
        ef = mo_surf.eltset_attribute("itetclr", 6)
        fs["front"] = ef.create_faceset(f"{base_name}_front.avs")
        return fs

    def createpts(
        self,
        crd: Literal["xyz", "rtz", "rtp"],
        npts: Tuple[int, int, int],
        mins: Tuple[float, float, float],
        maxs: Tuple[float, float, float],
        vc_switch=(1, 1, 1),
        rz_switch=(1, 1, 1),
        rz_value=(1, 1, 1),
        connect=False,
    ):
        """
        Create and Connect Points

        Parameters
        ----------
        crd : str
            Coordinate type of either 'xyz' (cartesian coordinates), 'rtz' (cylindrical coordinates), or 'rtp' (spherical coordinates).

        npts : Tuple[int, int, int]
            The number of points to create in line

        mins : Tuple[float, float, float]
            The starting value for each dimension.

        maxs : Tuple[float, float, float]
            The ending value for each dimension.

        vc_switch : Tuple[int, int, int]
            Determines if nodes represent vertices (1) or cell centers (0).

        rz_switch : Tuple[int, int, int]
            Determines if ratio zoning is used (1) or not (0).

        connect : bool
            Whether or not to connect points

        """

        cmd = "/".join(
            [
                "createpts",
                crd,
                ",".join([str(v) for v in npts]),
                ",".join([str(v) for v in mins]),
                ",".join([str(v) for v in maxs]),
                ",".join([str(v) for v in vc_switch]),
                ",".join([str(v) for v in rz_switch]),
                ",".join([str(v) for v in rz_value]),
            ]
        )
        self.sendcmd(cmd)

        if connect:
            if self.elem_type[0].startswith(("tri", "tet")):
                cmd = "/".join(["connect", "noadd"])
            else:
                cmd = "/".join(
                    [
                        "createpts",
                        "brick",
                        crd,
                        ",".join([str(v) for v in npts]),
                        "1,0,0",
                        "connect",
                    ]
                )
            self.sendcmd(cmd)

    def createpts_dxyz(
        self,
        dxyz: Tuple[float, float, float],
        mins: Tuple[float, float, float],
        maxs: Tuple[float, float, float],
        clip: Literal["under", "over"]
        | Tuple[
            Literal["under", "over"], Literal["under", "over"], Literal["under", "over"]
        ] = "under",
        hard_bound: Literal["min", "max"]
        | Tuple[
            Literal["min", "max"], Literal["min", "max"], Literal["min", "max"]
        ] = "min",
        rz_switch=(1, 1, 1),
        rz_value=(1, 1, 1),
        connect=True,
    ):
        """
        Create and Connect Points to create an orthogonal hexahedral mesh. The
        vertex spacing is based on dxyz and the mins and maxs specified. mins
        (default, see hard_bound option) or maxs will be adhered to, while maxs
        (default) or mins will be modified based on the clip option to be
        truncated at the nearest value 'under' (default) or 'over' the range
        maxs-mins. clip and hard_bound options can be mixed by specifying tuples
        (see description below).

        Parameters
        ----------
        dxyz : Tuple[float, float, float]
            The spacing between points in x, y, and z directions.

        mins : Tuple[float, float, float]
            The starting value for each dimension.

        maxs : Tuple[float, float, float]
            The ending value for each dimension.

        clip : str or Tuple[str, str, str], optional
            How to handle bounds if range does not divide by dxyz, either 'under' or 'over' range, by default 'under'

        hard_bound : str or Tuple[str, str, str], optional
            Whether to use the "min" or "max" as the hard constraint on dimension

        rz_switch : Tuple[int, int, int], optional
            Determines if ratio zoning is used (1) or not (0), by default (1, 1, 1)

        rzvalue : Tuple[int, int, int], optional
            Ratio zoning values, by default (1, 1, 1)

        connect : bool, optional
            Whether or not to connect points, by default True

        """
        if isinstance(hard_bound, str):
            hard_bound = (hard_bound, hard_bound, hard_bound)
        if isinstance(clip, str):
            clips = [clip, clip, clip]
        else:
            clips = clip
        dxyz_ar = np.array(dxyz)
        mins_ar = np.array(mins)
        maxs_ar = np.array(maxs)
        dxyz_ar[dxyz_ar == 0] = 1
        npts = np.zeros_like(dxyz_ar).astype("int")
        for i, cl in enumerate(clips):
            if cl == "under":
                npts[i] = int(np.floor((maxs_ar[i] - mins_ar[i]) / dxyz_ar[i]))
            elif cl == "over":
                npts[i] = int(np.ceil((maxs_ar[i] - mins_ar[i]) / dxyz_ar[i]))
        for i, bnd in enumerate(hard_bound):
            if bnd == "min":
                maxs_ar[i] = mins_ar[i] + npts[i] * dxyz_ar[i]
            elif bnd == "max":
                mins_ar[i] = maxs_ar[i] - npts[i] * dxyz_ar[i]
        npts += 1
        vc_switch = (1, 1, 1)  # always vertex nodes for dxyz method
        self.createpts(
            "xyz",
            npts.tolist(),
            mins_ar.tolist(),
            maxs_ar.tolist(),
            vc_switch,
            rz_switch,
            rz_value,
            connect=connect,
        )

    def createpts_line(
        self,
        npts: int,
        mins: Tuple[float, float, float],
        maxs: Tuple[float, float, float],
        vc_switch=(1, 1, 1),
        rz_switch=(1, 1, 1),
    ):
        """
        Create and Connect Points in a line

        Parameters
        ----------
        npts : int
            The number of points to create in line

        mins : Tuple[float, float, float]
            The starting value for each dimension.

        maxs : Tuple[float, float, float]
            The ending value for each dimension.

        vc_switch : Tuple[int, int, int], optional
            Determines if nodes represent vertices (1) or cell centers (0), by default (1, 1, 1)

        rz_switch : Tuple[int, int, int], optional
            Determines true or false (1 or 0) for using ratio zoning values, by default (1, 1, 1)

        """

        cmd = "/".join(
            [
                "createpts",
                "line",
                str(npts),
                " ",
                " ",
                ",".join([str(v) for v in mins + maxs]),
                ",".join([str(v) for v in vc_switch]),
                ",".join([str(v) for v in rz_switch]),
            ]
        )
        self.sendcmd(cmd)

    def createpts_brick(
        self,
        crd: Literal["xyz", "rtz", "rtp"],
        npts: Tuple[int, int, int],
        mins: Tuple[float, float, float],
        maxs: Tuple[float, float, float],
        vc_switch=(1, 1, 1),
        rz_switch=(1, 1, 1),
        rz_vls=(1, 1, 1),
    ):
        """
        Create and Connect Points

        Creates a grid of points in the mesh object and connects them.

        Parameters
        ----------
        crd : str
            Coordinate type of either 'xyz' (cartesian coordinates), 'rtz' (cylindrical coordinates), or 'rtp' (spherical coordinates).

        npts : Tuple[int, int, int]
            The number of points to create in each dimension.

        mins : Tuple[float, float, float]
            The starting value for each dimension.

        maxs : Tuple[float, float, float]
            The ending value for each dimension.

        vc_switch : Tuple[int, int, int]
            Determines if nodes represent vertices (1) or cell centers (0).

        rz_switch : Tuple[int, int, int]
            Determines if ratio zoning is used (1) or not (0).

        rz_vls : Tuple[int, int, int]
            Ratio zoning values.

        """

        ni, nj, nk = map(str, npts)
        xmn, ymn, zmn = map(str, mins)
        xmx, ymx, zmx = map(str, maxs)
        iiz, ijz, ikz = map(str, vc_switch)
        iirat, ijrat, ikrat = map(str, rz_switch)
        xrz, yrz, zrz = map(str, rz_vls)

        t = (crd, ni, nj, nk, xmn, ymn, zmn, xmx, ymx, zmx)
        t = t + (iiz, ijz, ikz, iirat, ijrat, ikrat, xrz, yrz, zrz)
        cmd = (
            "createpts/brick/%s/%s,%s,%s/%s,%s,%s/%s,%s,%s/%s,%s,%s/"
            + "%s,%s,%s/%s,%s,%s"
        )
        self.sendcmd(cmd % t)

    def createpts_median(self):
        self.sendcmd("createpts/median")

    def subset(
        self,
        mins: Tuple[float, float, float],
        maxs: Tuple[float, float, float],
        geom: Literal["xyz", "rtz", "rtp"] = "xyz",
    ):
        """
        Return Mesh Object Subset

        Creates a new mesh object that contains only a geometric subset defined
        by mins and maxs.

        Parameters
        ----------
        mins : Tuple[float, float, float]
            Coordinate of one of the shape's defining points.
                xyz (Cartesian):   (x1, y1, z1)
                rtz (Cylindrical): (radius1, theta1, z1)
                rtp (Spherical):   (radius1, theta1, phi1)

        maxs : Tuple[float, float, float]
            Coordinate of one of the shape's defining points.
                xyz (Cartesian):   (x2, y2, z2)
                rtz (Cylindrical): (radius2, theta2, z2)
                rtp (Spherical):   (radius2, theta2, phi2)

        geom : str, optional
            Type of geometric shape: 'xyz' (spherical), 'rtz' (cylindrical), 'rtp' (spherical)

        Returns
        -------
        MO

        Examples
        --------
        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()

        Create a mesh object.
        >>> mo = lg.create()
        >>> mo.createpts_brick("xyz", (5, 5, 5), (0, 0, 0), (5, 5, 5))

        Take the subset from (3,3,3)
        >>> subset = mo.subset((3, 3, 3), (5, 5, 5))

        >>> lg.close()

        """

        lg = self.lg
        new_mo = lg.copy(self)
        sub_pts = new_mo.pset_geom(mins, maxs, geom=geom)
        rm_pts = new_mo.pset_bool([sub_pts], boolean="not")

        new_mo.rmpoint_pset(rm_pts)
        return new_mo

    def quadxy(
        self,
        nnodes: Tuple[int, int, int],
        pts: Iterable[Tuple[float, float, float]],
        connect=True,
    ):
        """
        Define an arbitrary, logical quad of points in 3D space
        with nnodes(x,y,z) nodes. By default, the nodes will be connected.

        Parameters
        ----------
        nnode : Tuple[int, int, int]
            The number of nodes to create in each dimension.
                One value must == 1 and the other two must be > 1.

        pts : Iterable[Tuple[float, float, float]]
            The four corners of the quad surface, defined in counter
            clockwise order (the normal to the quad points is defined
            using the right hand rule and the order of the points).

        connect : bool, optional
            Connect the points, by default

        Examples
        --------
        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()

        Create a mesh object.
        >>> qua = lg.create("qua")

        Define 4 points in correct order
        >>> p1 = (0.0, 200.0, -400.0)
        >>> p2 = (0.0, -200.0, -400.0)
        >>> p3 = (140.0, -200.0, 0.0)
        >>> p4 = (118.0, 200.0, 0.0)
        >>> pts = [p1, p2, p3, p4]

        Define nnodes
        >>> nnodes = (29, 1, 82)

        Create and connect skewed plane
        >>> qua.quadxy(nnodes, pts)

        >>> lg.close()

        """
        self.select()
        quadpts = [n for n in nnodes if n != 1]
        assert len(quadpts) == 2, "nnodes must have one value == 1 and two values > 1"  # noqa: S101

        cmd = ["quadxy"]
        cmd.extend(map(str, quadpts))

        for v in pts:
            assert len(v) == 3, "vectors must be of length 3 (x,y,z)"  # noqa: S101
            cmd.append(",".join(list(map(str, v))))

        self.sendcmd("/".join(cmd))

        if connect:
            cmd = [
                "createpts",
                "brick",
                "xyz",
                ",".join(map(str, nnodes)),
                "1,0,0",
                "connect",
            ]

            self.sendcmd("/".join(cmd))

    def quadxyz(
        self,
        nnodes: Tuple[int, int, int],
        pts: List[Tuple[float, float, float]],
        connect=True,
    ):
        """
        Define an arbitrary and logical set of points in 3D (xyz) space.
        The set of points will be connected into hexahedrons by default. Set 'connect=False' to prevent connection.

        Parameters
        ----------
        nnodes : Tuple[int, int, int]
            The number of nodes to create in each dimension.
            The number of points will be 1 more than the number of elements in each dimension.

        pts : List[Tuple[float, float, float]]
            The eight corners of the hexahedron.
            The four bottom corners are listed first, then the four top corners.
            Each set of corners (bottom and top) are defined in counter-clockwise
            order (the normal to the quad points is defined using the right hand rule and the order of the points).

        connect : bool, optional
            Connect the points, by default

        Examples
        --------
        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()

        Create a mesh object.
        >>> hex = lg.create()

        Define 4 bottom points in correct order
        >>> p1 = (0.0, 0.0, 0.0)
        >>> p2 = (1.0, 0.0, 0.02)
        >>> p3 = (1.0, 1.0, 0.0)
        >>> p4 = (0.0, 1.0, 0.1)

        Define 4 top points in correct order
        >>> p5 = (0.0, 0.0, 1.0)
        >>> p6 = (1.0, 0.0, 1.0)
        >>> p7 = (1.0, 1.0, 1.0)
        >>> p8 = (0.0, 1.0, 1.1)
        >>> pts = [p1, p2, p3, p4, p5, p6, p7, p8]

        Define nnodes
        >>> nnodes = (3, 3, 3)

        Create and connect skewed hex mesh
        >>> hex.quadxyz(nnodes, pts)

        Dump mesh
        >>> hex.dump("quadxyz_test.gmv")

        >>> lg.close()

        """
        self.select()
        if len(pts) != 8:
            raise ValueError("pts must contain eight sets of points")
        cmd = ["quadxyz", ",".join(map(str, nnodes))]
        for v in pts:
            cmd.append(",".join(list(map(str, v))))

        self.sendcmd("/".join(cmd))

        if connect:
            cmd = [
                "createpts",
                "brick",
                "xyz",
                ",".join(map(str, nnodes)),
                "1,0,0",
                "connect",
            ]

            self.sendcmd("/".join(cmd))

    def rzbrick(
        self,
        n_ijk: Tuple[int, int, int],
        connect=True,
        stride=(1, 0, 0),
        coordinate_space: Literal["xyz", "rtz", "rtp"] = "xyz",
    ):
        """
        Builds a brick mesh and generates a nearest neighbor connectivity matrix

        Currently only configured for this flavor of syntax:

            rzbrick/xyz|rtz|rtp/ni,nj,nk/pset,get,name/connect/

        Use this option with quadxyz to connect logically rectangular grids.

        Parameters
        ----------
        n_ij : Tuple[int, int, int]
            The number of points to create in each direction.

        connect : bool, optional
            Connect the points, by default

        stride : Tuple[int, int, int], optional
            Stride to select, by default (1, 0, 0)

        coordinate_space : str, optional
            The coordinate space to use, by default "xyz"

        """
        self.select()
        cmd = f"rzbrick/{coordinate_space}"

        for v in [n_ijk, stride]:
            cmd += "/" + ",".join(list(map(str, v)))

        if connect:
            cmd += "/connect"

        self.sendcmd(cmd)

    def grid2grid(
        self,
        ioption: Literal[
            "quadtotri2",
            "prismtotet3",
            "quattotri4",
            "pyrtotet4",
            "hextotet5",
            "hextotet6",
            "prismtotet14",
            "prismtotet18",
            "hextotet24",
            "tree_to_fe",
        ],
        name: Optional[str] = None,
    ):
        """
        Convert a mesh with one element type to a mesh with another

        Parameters
        ----------
        ioption : str
            Type of conversion:
                quadtotri2      quad to 2 triangles, no new points.
                prismtotet3     prism to 3 tets, no new points.
                quadtotri4      quad to 4 triangles, with one new point.
                pyrtotet4       pyramid to 4 tets, with one new point.
                hextotet5       hex to 5 tets, no new points.
                hextotet6       hex to 6 tets, no new points.
                prismtotet14    prism to 14 tets, four new points (1 + 3 faces).
                prismtotet18    prism to 18 tets, six new points (1 + 5 faces).
                hextotet24      hex to 24 tets, seven new points (1 + 6 faces).
                tree_to_fe      quadtree or octree grid to grid with no parent-type elements.

        name : str, optional
            Internal Lagrit name of new mesh object.

        Returns
        -------
        MO

        """
        if name is None:
            name = new_name("mo", self.lg.core.mo_names())
        cmd = "/".join(["grid2grid", ioption, name, self.name])
        self.sendcmd(cmd)
        return MO(name, self.lg)

    def connect(
        self,
        algorithm: Optional[Literal["delaunay", "noadd", "check_interface"]] = None,
        option: Optional[str] = None,
        stride: Optional[Tuple[int, int, int]] = None,
        big_tet_coords: Iterable[Tuple[float, float, float]] = [],  # noqa: B006
    ):
        """
        Connect the nodes into a Delaunay tetrahedral or triangle grid.

        Parameters
        ----------
        algorithm : str, optional
            Type of connect: delaunay, noadd, or check_interface

        option : str, optional
            type of connect: noadd, or check_interface

        stride : Tuple[int, int, int], optional
            tuple of (first, last, stride) of points

        """
        cmd = ["connect"]
        if algorithm is not None:
            cmd.append(algorithm)
        if stride is not None and algorithm == "delaunay":
            cmd += [",".join(map(str, stride))]
            for b in big_tet_coords:
                cmd += [",".join(map(str, b))]
        if option is not None:
            cmd += [option]
        cmd = "/".join(cmd)
        self.sendcmd(cmd)

    def connect_delaunay(
        self,
        option: Optional[str] = None,
        stride: Optional[Tuple[int, int, int]] = None,
        big_tet_coords: Iterable[Tuple[float, float, float]] = [],  # noqa: B006
    ):
        """
        Connect the nodes into a Delaunay tetrahedral or triangle grid without adding nodes.
        """
        mo_tmp = self.copypts()
        mo_tmp.setatt("imt", 1)
        mo_tmp.setatt("itp", 0)
        mo_tmp.rmpoint_compress(filter_bool=True)
        mo_tmp.connect(
            algorithm="delaunay",
            option=option,
            stride=stride,
            big_tet_coords=big_tet_coords,
        )
        self.sendcmd("/".join(["cmo", "move", self.name, mo_tmp.name]))

    def connect_noadd(self):
        """
        Connect the nodes into a Delaunay tetrahedral or triangle grid without adding nodes.
        """
        self.connect(algorithm="noadd")

    def connect_check_interface(self):
        """
        Connect the nodes into a Delaunay tetrahedral or triangle grid
        exhaustively checking that no edges of the mesh cross a material
        boundary.
        """
        self.connect(algorithm="check_interface")

    def copypts(
        self,
        elem_type: Literal[
            "tet", "hex", "pri", "pyr", "tri", "qua", "hyb", "lin", "triplane"
        ] = "tet",
        name: Optional[str] = None,
    ):
        """
        Copy points from mesh object to new mesh object

        Parameters
        ----------
        name : str, optional
            Internal Lagrit name of new mesh object.

        mesh_type : str, optional
            Mesh type for new mesh

        Returns
        -------
        MO

        """
        if name is None:
            name = new_name("mo", self.lg.core.mo_names())
        mo_new = self.lg.create(elem_type=elem_type, name=name)
        self.sendcmd("/".join(["copypts", mo_new.name, self.name]))
        return mo_new

    def extrude(
        self,
        offset: float,
        offset_type="const",
        return_type="volume",
        direction: Optional[Tuple[float, float, float]] = None,
        name: Optional[str] = None,
    ) -> "MO":
        """
        Extrude mesh object to new mesh object
        This command takes the current mesh object (topologically 1d or 2d mesh (a line, a set of line
        segments, or a planar or non-planar surface)) and extrudes it into three
        dimensions along either the normal to the curve or surface (default),
        along a user defined vector, or to a set of points that the user has specified.
        If the extrusion was along the normal of the surface or along a user
        defined vector, the command can optionally find the external surface of
        the volume created and return that to the user.
        Refer to http://lagrit.lanl.gov/docs/commands/extrude.html for more details on arguments.

        Parameters
        ----------
        offfset : float
            Distance to extrude

        offset_type : str, optional
            either const or min (interp will be handled in the PSET class in the future)

        return_type : str, optional
            either volume for entire mesh or bubble for just the external surface

        direction : Tuple[float, float, float], optional
            Direction to extrude in, defaults to normal of the object

        name : str, optional
            Name to use within lagrit for the created mesh object

        Returns
        -------
        MO

        """
        if name is None:
            name = new_name("mo", self.lg.core.mo_names())
        cmd = ["extrude", name, self.name, offset_type, str(offset), return_type]
        if direction is not None:
            cmd.append(",".join(map(str, direction)))
        self.sendcmd("/".join(cmd))
        return MO(name, self.lg)

    def refine_to_object(
        self,
        mo: "MO",
        level: int = 1,
        imt: Optional[int] = None,
        prd_choice: Optional[int] = None,
    ):
        """
        Refine mesh at locations that intersect another mesh object

        Parameters
        ----------
        mo : MO
            Mesh object to intersect with current mesh object to determine where to refine

        level : int, optional
            max level of refinement

        imt : int, optional
            Value to assign to imt (LaGriT material type attribute)

        prd_choice : int
            directions of refinement

        """
        for _ in range(level):
            attr_name = self.intersect_elements(mo)

            if level > 1:
                e_attr = self.eltset_attribute(attr_name, 0, boolstr="gt")
                e_level = self.eltset_attribute("itetlev", level, boolstr="lt")
                e_refine = self.eltset_bool([e_attr, e_level], boolstr="inter")
                e_attr.delete()
                e_level.delete()
            else:
                e_refine = self.eltset_attribute(attr_name, 0, boolstr="gt")

            if prd_choice is not None:
                p_refine = e_refine.pset()
                p_refine.refine(prd_choice=prd_choice)
                p_refine.delete()
            else:
                e_refine.refine()

            mo.rmpoint_eltset(e_refine)
            e_refine.delete()

        if imt is not None:
            attr_name = self.intersect_elements(mo)
            e_attr = self.eltset_attribute(attr_name, 0, boolstr="gt")
            p = e_attr.pset()
            p.setatt("imt", imt)
            p.delete()

    def intersect_elements(self, mo: "MO", attr_name: Optional[str] = None):
        """
        This command takes two meshes and creates an element-based attribute in mesh1
        that contains the number of elements in mesh2 that intersected the respective
        element in mesh1. We define intersection as two elements sharing any common point.

        Parameters
        ----------
        mo : MO
            Mesh object to intersect with current mesh object to determine where to refine

        attr_name : str, optional
            Name of attribute to create

        Returns
        -------
        str of attribute name

        """
        attr_name = attr_name if attr_name else "attr00"
        self.sendcmd("/".join(["intersect_elements", self.name, mo.name, attr_name]))
        return attr_name

    def extract_surfmesh(
        self,
        name: Optional[str] = None,
        stride=(1, 0, 0),
        reorder=False,
        resetpts_itp=True,
        external=False,
    ):
        return self.lg.extract_surfmesh(
            name=name,
            cmo_in=self,
            stride=stride,
            reorder=reorder,
            resetpts_itp=resetpts_itp,
            external=external,
        )

    def interpolate(
        self,
        method: Literal["map", "voronoi", "continuous", "default"],
        attsink: str,
        cmosrc: "MO",
        attsrc: str,
        stride=(1, 0, 0),
        tie_option: Optional[Literal["tiemin", "tiemax"]] = None,
        flag_option: Optional[Literal["plus1", "nearest"] | int | float] = None,
        keep_option: Optional[Literal["keepatt", "delatt"]] = None,
        interp_function: Optional[
            Literal[
                "linear",
                "asinh",
                "log",
                "copy",
                "sequence",
                "min",
                "incmin",
                "max",
                "incmax",
                "and",
                "or",
            ]
            | int
            | float
        ] = None,
    ):
        """
        Interpolate values from attribute attsrc from mesh object cmosrc to current mesh object
        """
        stride = [str(v) for v in stride]
        cmd = [
            "interpolate",
            method,
            self.name,
            attsink,
            ",".join(stride),
            cmosrc.name,
            attsrc,
        ]
        if tie_option is not None:
            cmd.append(tie_option)
        if flag_option is not None:
            cmd.append(str(flag_option))
        if keep_option is not None:
            cmd.append(keep_option)
        if interp_function is not None:
            cmd.append(str(interp_function))
        self.sendcmd("/".join(cmd))

    def interpolate_continuous(
        self,
        attsink: str,
        cmosrc: "MO",
        attsrc: str,
        stride=(1, 0, 0),
        interp_function: Optional[str] = None,
        nearest: Optional[str] = None,
    ):
        stride = [str(v) for v in stride]
        cmd = [
            "intrp",
            "continuous",
            self.name + " " + attsink,
            ",".join(stride),
            cmosrc.name + " " + attsrc,
        ]
        if nearest is not None:
            cmd += ["nearest", nearest]
        if interp_function is not None:
            cmd.append(interp_function)
        print("/".join(cmd))
        self.sendcmd("/".join(cmd))

    def copy(self, name: Optional[str] = None):
        """
        Copy mesh object
        """
        if name is None:
            name = new_name("mo", self.lg.core.mo_names())
        self.sendcmd("/".join(["cmo/copy", name, self.name]))
        return MO(name, self.lg)

    def stack_layers(
        self,
        filelist: List[str],
        file_type="avs",
        nlayers: Optional[List[int]] = None,  # ref_num
        matids: Optional[List[int]] = None,
        xy_subset: Optional[
            Tuple[float, float, float, float]
        ] = None,  # minx, miny, maxx, maxy
        buffer_opt: Optional[float] = None,
        truncate_opt: Optional[int] = None,
        pinchout_opt: Optional[float] = None,
        dpinchout_opt: Optional[Tuple[float, float]] = None,
        flip_opt=False,
    ):
        if nlayers is None:
            refine_num = [""] * (len(filelist) - 1)
        else:
            refine_num = [str(n) for n in nlayers]
        if matids is None:
            mat_num = ["1"] * len(filelist)
        else:
            mat_num = [str(m) for m in matids]
        cmd = ["stack", "layers", file_type]
        if xy_subset is not None:
            cmd.append(",".join([str(v) for v in xy_subset]))

        cmd.append(f"{filelist[0]} {mat_num[0]}")
        for file, matid, n_refine in zip(filelist[1:], mat_num[1:], refine_num):
            cmd.append(f"{file} {matid} {n_refine}")
        cmd.append(f"{filelist[-1]} {mat_num[-1]} {refine_num[-1]}")

        if flip_opt is True:
            cmd.append("flip")
        if buffer_opt is not None:
            cmd.append(f"buffer {buffer_opt}")
        if truncate_opt is not None:
            cmd.append(f"trunc {truncate_opt}")
        if pinchout_opt is not None:
            cmd.append(f"pinch {pinchout_opt}")
        if dpinchout_opt is not None:
            cmd.append(f"dpinch {dpinchout_opt[0]}")
            cmd.append(f"dmin {dpinchout_opt[1]}")

        self.lg.sendcmd("/".join(cmd))

    def stack_fill(self, name: Optional[str] = None):
        if name is None:
            name = new_name("mo", self.lg.core.mo_names())
        self.sendcmd("/".join(["stack/fill", name, self.name]))
        return MO(name, self.lg)

    def math(
        self,
        operation: str,
        attsink: str,
        value: Optional[float] = None,
        stride: Tuple[int, int, int] | Tuple[str, str, str] = (1, 0, 0),
        cmosrc: Optional["MO"] = None,
        attsrc: Optional[str] = None,
    ):
        if cmosrc is None:
            cmosrc = self
        if attsrc is None:
            attsrc = attsink
        cmd = [
            "math",
            operation,
            self.name,
            attsink,
            ",".join([str(v) for v in stride]),
            cmosrc.name,
            attsrc,
        ]
        if value is not None:
            cmd += [str(value)]
        self.sendcmd("/".join(cmd))

    def settets(
        self,
        method: Optional[
            Literal[
                "color_tets",
                "parents",
                "geometry",
                "color_points",
                "newtets",
                "normal",
            ]
        ] = None,
    ):
        if method is None:
            self.sendcmd("settets")
        else:
            self.sendcmd("settets/" + method)

    def triangulate(
        self, order: Literal["clockwise", "counterclockwise"] = "clockwise"
    ):
        """
        triangulate will take an ordered set of nodes in the current 2d mesh object that define a perimeter of a polygon and create a trangulation of the polygon.  The nodes are assumed to lie in the xy plane; the z coordinate is ignored.  No checks are performed to verify that the nodes define a legal perimeter (i.e. that segments of the perimeter do not cross).  The code will connect the last node to the first node to complete the perimeter.

        This code support triangulation of self-intersecting polygons (polygon with holes), assuming that the order of the nodes are correct. Moreover the connectivity of the polyline must also be defined correctly. No checks are made.

        One disadvantage of the algorithm for triangulating self-intersecting polygons is that it does not always work. For example, if the holes have complicated shapes, with many concave vertices, the code might fail. In this case, the user may try to rotate the order of the nodes:
            NODE_ID:
                1 -> 2
                2 -> 3
                ...
                N -> 1

        Parameters
        ----------
        order : str, optional
            direction of point ordering, by default "clockwise"

        Examples
        --------
        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()

        Define polygon points in clockwise direction and create tri mesh object
        >>> coords = [
        ...     (0.0, 0.0, 0.0),
        ...     (0.0, 1000.0, 0.0),
        ...     (2200.0, 200.0, 0.0),
        ...     (2200.0, 0.0, 0.0),
        ... ]
        >>> motri = lg.tri_mo_from_polyline(coords)

        Triangulate polygon
        >>> motri.triangulate()
        >>> motri.setatt("imt", 1)
        >>> motri.setatt("itetclr", 1)

        Refine mesh with successively smaller edge length constraints
        >>> edge_length = [1000, 500, 250, 125, 75, 40, 20, 15]
        >>> for i, l in enumerate(edge_length):
        ...     motri.resetpts_itp()
        ...     motri.refine(refine_option='rivara',refine_type='edge',values=[l],inclusive_flag='inclusive')
        ...     motri.smooth()
        ...     motri.recon(0)

        Provide additional smoothing after the last refine
        >>> for i in range(5):
        ...     motri.smooth()
        ...     motri.recon(0)

        Create delaunay mesh and clean up
        >>> motri.tri_mesh_output_prep()

        Dump fehm files
        >>> motri.dump("fehm", "nk_mesh00")

        >>> lg.close()

        """
        self.sendcmd("triangulate/" + order)

    def refine(
        self,
        refine_option="constant",
        field=" ",
        interpolation=" ",
        refine_type="element",
        stride=(1, 0, 0),
        values=[1.0],  # noqa: B006
        inclusive_flag="exclusive",
        prd_choice: Optional[int] = None,
    ):
        cmd = [
            "refine",
            refine_option,
            field,
            interpolation,
            refine_type,
            " ".join([str(v) for v in stride]),
            "/".join([str(v) for v in values]),
            inclusive_flag,
        ]
        if prd_choice is not None:
            cmd.append(f"amr {prd_choice}")
        self.lg.sendcmd("/".join(cmd))

    def regnpts(
        self,
        geom: Literal["xyz"] | Literal["rtz"] | Literal["rtp"],
        ray_points: Tuple[  # xyz
            Tuple[float, float, float],
            Tuple[float, float, float],
            Tuple[float, float, float],
        ]
        | Tuple[Tuple[float, float, float], Tuple[float, float, float]]  # rtz
        | Tuple[Tuple[float, float, float]],  # rpt
        region: "Region",
        ptdist: Literal["inside", "in", "out", "outside", "both"] | float | int,
        stride: Tuple[float, float, float] | "PSet" = (1, 0, 0),
        irratio: float = 0,
        rrz: float = 0,
        maxpenetr=False,
    ):
        """
        Generates points in a region previously defined by the region command.
        The points are generated by shooting rays through a user specified set of points from a
        plane and finding the intersection of each ray with the surfaces that define the region.

        Parameters
        ----------
        geom : str
            Type of geometric shape: 'xyz' (spherical), 'rtz' (cylindrical), 'rtp' (spherical)

        ray_points : Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]
            Three points that define plane which rays emante from

        region : Region
            Region to generate points within

        ptdist : str
            Parameter that determines point distribution pattern

        stride : Tuple[int, int, int], optional
            Points to shoot rays through

        irratio : float, optional
            Parameter that determines point distribution pattern

        rrz : float, optional
            Ratio zoning value

        maxpenetr : bool, optional
            Maximum distance along ray that points will be distributed

        Examples
        --------
        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()

        >>> p1 = (30.0, 0.0, 0.0)
        >>> p2 = (30.0, 1.0, 0.0)
        >>> p3 = (30.0, 1.0, 0.1)
        >>> pts = [p1, p2, p3]

        >>> npts = (3, 3, 3)
        >>> mins = (0, 0, 0)
        >>> maxs = (10, 10, 10)

        >>> mesh = lg.createpts("xyz", npts, mins, maxs, "hex", connect=False)
        >>> rayend = mesh.pset_geom(mins, maxs, ctr=(5, 5, 5), geom="xyz")
        >>> mesh.rmpoint_compress(filter_bool=True)
        >>> eighth = mesh.surface_box(mins, (5, 5, 5))
        >>> boolstr2 = f"gt {eighth.name}"
        >>> reg2 = mesh.region(boolstr2)
        >>> mesh.regnpts("xyz", pts, reg2, 1000, stride=rayend)

        >>> lg.close()

        """
        cmd = ["regnpts", region.name, str(ptdist)]

        if isinstance(stride, PSet):
            cmd.append(f"pset get {stride.name}")
        else:
            cmd.append(",".join([str(v) for v in stride]))

        cmd.append(geom)

        if geom == "xyz":
            if len(ray_points) != 3:
                raise ValueError("ray_points must contain three sets of points")
            for p in ray_points:
                cmd.append(",".join(list(map(str, p))))
        elif geom == "rtz":
            if len(ray_points) != 2:
                raise ValueError("ray_points must contain two sets of points")
            for p in ray_points:
                cmd.append(",".join(list(map(str, p))))
        elif geom == "rtp":
            if len(ray_points) != 1:
                raise ValueError("ray_points must contain one set of points")
            for p in ray_points:
                cmd.append(",".join(list(map(str, p))))

        cmd.append(str(irratio))
        cmd.append(str(rrz))
        if maxpenetr:
            cmd.append("maxpenetr")

        self.sendcmd("/".join(cmd))

    def setpts(self, no_interface=False, closed_surfaces=False):
        """
        Set point types and imt material by calling surfset and regset routines.
        """

        cmd = ["setpts"]
        if no_interface and closed_surfaces:
            raise ValueError("no_interface and closed_surfaces are mutually exclusive")

        if no_interface:
            cmd.append("no_interface")
        elif closed_surfaces:
            cmd.extend(["closed_surfaces", "reflect"])

        self.sendcmd("/".join(cmd))

    def smooth(self, *args, **kwargs):
        if "algorithm" not in kwargs:
            self.sendcmd("smooth")
        else:
            cmd = ["smooth", "position", kwargs["algorithm"]]
            for a in args:
                cmd.append(a)
            self.sendcmd("/".join(cmd))

    def recon(self, option: Literal[0, 1] = 0, damage="", checkaxy=False):
        cmd = ["recon", str(option), str(damage)]
        if checkaxy:
            cmd.append("checkaxy")
        self.sendcmd("/".join(cmd))

    def filter(
        self,
        stride=(1, 0, 0),
        tolerance: Optional[float] = None,
        boolean: Optional[Literal["min", "max"]] = None,
        attribute: Optional[str] = None,
    ):
        stride = [str(v) for v in stride]
        cmd = ["filter", " ".join(stride)]
        if tolerance is not None:
            cmd.append(str(tolerance))
        if boolean is not None and attribute is not None:
            cmd.append(boolean)
            cmd.append(attribute)
        elif (boolean is None and attribute is not None) or (
            boolean is not None and attribute is None
        ):
            raise ValueError(
                "Error: Both boolean and attribute must be specified together"
            )

        self.sendcmd("/".join(cmd))

    def tri_mesh_output_prep(self):
        """
        Prepare tri mesh for output, remove dudded points,
        ensure delaunay volumes, etc.
        Combination of lagrit commands:
        filter/1 0 0
        rmpoint/compress
        recon/1
        resetpts/itp
        """
        self.filter()
        self.rmpoint_compress()
        self.recon(1)
        self.resetpts_itp()

    def surface(self, name: Optional[str] = None, ibtype="reflect"):
        if name is None:
            name = new_name("s", self.surfaces.keys())
        cmd = "/".join(["surface", name, ibtype, "sheet", self.name])
        self.sendcmd(cmd)
        self.surfaces[name] = Surface(name, self)
        return self.surfaces[name]

    def surface_box(
        self,
        mins: Tuple[float, float, float],
        maxs: Tuple[float, float, float],
        name: Optional[str] = None,
        ibtype="reflect",
    ):
        if name is None:
            name = new_name("s", self.surfaces.keys())
        cmd = "/".join(
            [
                "surface",
                name,
                ibtype,
                "box",
                ",".join([str(v) for v in mins]),
                ",".join([str(v) for v in maxs]),
            ]
        )
        self.sendcmd(cmd)
        self.surfaces[name] = Surface(name, self)
        return self.surfaces[name]

    def surface_cylinder(
        self,
        coord1: Tuple[float, float, float],
        coord2: Tuple[float, float, float],
        radius: float,
        name: Optional[str] = None,
        ibtype="reflect",
    ):
        if name is None:
            name = new_name("s", self.surfaces.keys())
        cmd = "/".join(
            [
                "surface",
                name,
                ibtype,
                "cylinder",
                ",".join([str(v) for v in coord1]),
                ",".join([str(v) for v in coord2]),
                str(radius),
            ]
        )
        self.sendcmd(cmd)
        self.surfaces[name] = Surface(name, self)
        return self.surfaces[name]

    def surface_plane(
        self,
        coord1: Tuple[float, float, float],
        coord2: Tuple[float, float, float],
        coord3: Tuple[float, float, float],
        name: Optional[str] = None,
        ibtype="reflect",
    ):
        if name is None:
            name = new_name("s", self.surfaces.keys())
        cmd = [
            "surface",
            name,
            ibtype,
            "plane",
            ",".join(map(str, coord1)),
            ",".join(map(str, coord2)),
            ",".join(map(str, coord3)),
        ]
        self.sendcmd("/".join(cmd))
        self.surfaces[name] = Surface(name, self)
        return self.surfaces[name]

    def region(
        self,
        boolstr: str,
        name: Optional[str] = None,
    ):
        """
        Create region using boolean string

        Paramaters
        ----------
        boolstr : str
            String of boolean operations

        name : str, optional
            Internal lagrit name for mesh object

        Returns
        -------
        Region

        Examples
        --------
        >>> import pylagrit
        >>> lg = pylagrit.PyLaGriT()

        >>> mesh = lg.create()
        >>> mins = (0, 0, 0)
        >>> maxs = (5, 5, 5)
        >>> eighth = mesh.surface_box(mins, maxs)
        >>> boolstr1 = f"le {eighth.name}"
        >>> boolstr2 = f"gt {eighth.name}"
        >>> reg1 = mesh.region(boolstr1)
        >>> reg2 = mesh.region(boolstr2)
        >>> mreg1 = mesh.mregion(boolstr1)
        >>> mreg2 = mesh.mregion(boolstr2)
        >>> mesh.createpts_brick("xyz", (10, 10, 10), (0, 0, 0), (10, 10, 10))
        >>> mesh.rmregion(reg1)

        >>> lg.close()

        """
        if name is None:
            name = new_name("r", self.regions.keys())
        cmd = "/".join(["region", name, boolstr])
        self.sendcmd(cmd)
        self.regions[name] = Region(name, self)
        return self.regions[name]

    def mregion(self, boolstr: str, name: Optional[str] = None):
        """
        Create mregion using boolean string

        Paramaters
        ----------
        boolstr : str
            String of boolean operations

        name : str, optional
            Internal lagrit name for mesh object

        Returns
        -------
        MRegion

        """
        if name is None:
            name = new_name("mr", self.mregions.keys())
        cmd = "/".join(["mregion", name, boolstr])
        self.sendcmd(cmd)
        self.mregions[name] = MRegion(name, self)
        return self.mregions[name]

    def rmregion(
        self, region: "Region", rmpoints=True, filter_bool=False, resetpts_itp=True
    ):
        """
        Remove points that lie inside region

        Parameters
        ----------
        region : Region
            Region to remove points from

        """
        cmd = "/".join(["rmregion", region.name])
        self.sendcmd(cmd)
        if rmpoints:
            self.rmpoint_compress(filter_bool=filter_bool, resetpts_itp=resetpts_itp)

    def quality(
        self,
        *args,
        quality_type: Optional[
            Literal["aspect", "edge_ratio", "edge_min", "edge_max", "angle", "pcc"]
        ] = None,
        save_att=False,
    ):
        cmd = ["quality"]
        if quality_type is not None:
            cmd.append(quality_type)
            if save_att:
                cmd.append("y")
            for a in args:
                cmd.append(a)
        self.sendcmd("/".join(cmd))

    def rmmat(
        self,
        material_number: int,
        option: Literal["node", "element", "all"] = "all",
        exclusive=False,
    ):
        """
        This routine is used to remove points that are of a specified material value
        (itetclr for elements or imt for nodes). Elements with the specified material
        value are flagged by setting the element material type negative. They are not
        removed from the mesh object.

        Parameters
        ----------
        material_number : int
            Number of material

        option : str, optional
            'node' removes nodes with imt=material_number
            'element' removes elements with itetclr=material_number
            'all' removes nodes and elements with material_number equal to imt and itetclr, respectively

        exclusive : bool, optional
            if True, removes everything except nodes with imt=material and removes everything except elements with itetclr=material number.

        """
        cmd = ["rmmat", str(material_number), option]
        if exclusive:
            cmd.append("exclusive")
        self.sendcmd("/".join(cmd))


class Surface:
    """
    Surface
    """

    def __init__(self, name: str, parent: MO):
        self.name = name
        self._parent = parent

    def __repr__(self):
        return self.name

    def release(self):
        cmd = "surface/" + self.name + "/release"
        self._parent.sendcmd(cmd)
        del self._parent.surfaces[self.name]


class Region:
    """
    Region
    """

    def __init__(self, name: str, parent: MO):
        self.name = name
        self._parent = parent

    def __repr__(self):
        return str(self.name)

    def release(self):
        cmd = "region/" + self.name + "/release"
        self._parent.sendcmd(cmd)
        del self._parent.regions[self.name]

    def eltset(self, name: Optional[str] = None):
        if name is None:
            name = new_name("elt", self._parent.lg.core.eltset_names())
        cmd = "/".join(["eltset", name, "region", str(self)])
        self._parent.sendcmd(cmd)
        return EltSet(name, self._parent)

    def dump_pflotran_region(self, filename: str):
        elts = self.eltset(name=self.name)
        elts.dump_pflotran_region(filename)
        elts.delete()


class MRegion:
    """
    Material Region
    """

    def __init__(self, name: str, parent: MO):
        self.name = name
        self._parent = parent

    def __repr__(self):
        return str(self.name)

    def release(self):
        cmd = "mregion/" + self.name + "/release"
        self._parent.sendcmd(cmd)
        del self._parent.mregions[self.name]

    def eltset(self, name: Optional[str] = None):
        if name is None:
            name = new_name("elt", self._parent.lg.core.eltset_names())
        cmd = "/".join(["eltset", name, "mregion", str(self)])
        self._parent.sendcmd(cmd)
        return EltSet(name, self._parent)

    def dump_pflotran_region(self, filename: str):
        elts = self.eltset(name=self.name)
        elts.dump_pflotran_region(filename)
        elts.delete()


class EltSet:
    """
    EltSet
    """

    def __init__(self, name: str, parent: MO):
        self.name = name
        self.faceset: Optional[FaceSet] = None
        self._parent = parent

    def __repr__(self):
        return str(self.name)

    def delete(self):
        cmd = "eltset/" + self.name + "/delete"
        self._parent.sendcmd(cmd)

    def create_faceset(self, filename: Optional[str] = None):
        if filename is None:
            filename = "faceset_" + self.name + ".avs"
        motmpnm = new_name("mo_tmp", self._parent.lg.core.mo_names())
        self._parent.lg.sendcmd("/".join(["cmo/copy", motmpnm, self._parent.name]))
        self._parent.lg.sendcmd("/".join(["cmo/DELATT", motmpnm, "itetclr0"]))
        self._parent.lg.sendcmd("/".join(["cmo/DELATT", motmpnm, "itetclr1"]))
        self._parent.lg.sendcmd("/".join(["cmo/DELATT", motmpnm, "facecol"]))
        self._parent.lg.sendcmd("/".join(["cmo/DELATT", motmpnm, "idface0"]))
        self._parent.lg.sendcmd("/".join(["cmo/DELATT", motmpnm, "idelem0"]))
        self._parent.lg.sendcmd("eltset / eall / itetclr / ge / 0")
        self._parent.lg.sendcmd("eltset/edel/not eall " + self.name)
        self._parent.lg.sendcmd("rmpoint / element / eltset get edel")
        self._parent.lg.sendcmd("rmpoint / compress")
        self._parent.lg.sendcmd("/".join(["dump / avs2", filename, motmpnm, "0 0 0 2"]))
        self._parent.lg.sendcmd("cmo / delete /" + motmpnm)
        self.faceset = FaceSet(filename, self)
        return self.faceset

    def minmax(self, attname: Optional[str] = None, stride=(1, 0, 0)):
        self._parent.printatt(
            attname=attname, stride=stride, eltset=self, ptype="minmax"
        )

    def list(self, attname: Optional[str] = None, stride=(1, 0, 0)):
        self._parent.printatt(attname=attname, stride=stride, eltset=self, ptype="list")

    def refine(self):
        """
        Refine elements in the element set

        Examples
        --------
        >>> import pylagrit
        >>> import numpy as np
        >>> import sys
        >>> lg = pylagrit.PyLaGriT()

        >>> df = 0.0005  # Fault half aperture
        >>> lr = 7  # Levels of refinement
        >>> nx = 4  # Number of base mesh blocks in x direction
        >>> nz = 20  # Number of base mesh blocks in z direction
        >>> d_base = df * 2 ** (lr + 1)  # Calculated dimension of base block
        >>> w = d_base * nx  # Calculated width of model
        >>> d = d_base * nz  # Calculated depth of model

        Create discrete fracture mesh
        >>> dxyz = np.array([d_base, d_base, 0.0])
        >>> mins = np.array([0.0, -d, 0.0])
        >>> maxs = np.array([w, 0, 0])
        >>> mqua = lg.createpts_dxyz(
        ...     dxyz, mins, maxs, "qua", hard_bound=("min", "max", "min"), connect=True
        ... )

        >>> for i in range(lr):
        ...     prefine = mqua.pset_geom("xyz", mins-0.1,(0.0001,0.1,0))
        ...     erefine = prefine.eltset()
        ...     erefine.refine()
        ...     prefine.delete()
        ...     erefine.delete()
        >>> mtri = mqua.copypts("triplane")
        >>> mtri.connect()

        >>> mtri.tri_mesh_output_prep()
        >>> mtri.reorder_nodes(cycle="xic yic zic")
        >>> pfault = mtri.pset_geom("xyz", mins - 0.1, (0.0001, 0.1, 0))
        >>> psource = mtri.pset_geom("xyz", mins - 0.1, mins + 0.0001)
        >>> mtri.setatt("imt", 1)
        >>> pfault.setatt("imt", 10)
        >>> psource.setatt("imt", 20)

        >>> lg.close()

        """
        cmd = "/".join(["refine", "eltset", f"eltset,get,{self}"])
        self._parent.sendcmd(cmd)

    def pset(self, name: Optional[str] = None):
        """
        Create a pset from the points in an element set

        Parameters
        ----------

        name : str, optional
            Name of point set to be used within LaGriT

        Returns
        -------
        PSet
        """
        if name is None:
            name = new_name("p", self._parent.lg.core.pset_names())
        cmd = "/".join(["pset", name, "eltset", self.name])
        self._parent.sendcmd(cmd)
        return PSet(name, self._parent)

    def setatt(self, attname: str, value: int | float | str):
        cmd = "/".join(
            [
                "cmo/setatt",
                self._parent.name,
                attname,
                "eltset, get," + self.name,
                str(value),
            ]
        )
        self._parent.sendcmd(cmd)

    def dump_pflotran_region(self, filename: str):
        att_name = new_name("mask", self._parent.information()["attributes"].keys())
        self._parent.addatt(att_name, vtype="vint", length="nelements", value=0)
        self.setatt(att_name, 1)

        mask = cast(List[int], self._parent.obj.attr(att_name))

        if filename.endswith(".h5"):
            from h5py import File

            h5file = File(filename, mode="a")
            h5file.create_dataset(
                f"Regions/{self.name}/Cell Ids",
                data=np.arange(1, len(mask) + 1)[mask == 1],
            )
            h5file.close()
        else:
            with open(filename, "w") as f:
                for i, m in enumerate(mask):
                    if m == 1:
                        f.write(f"{i+1}\n")

        self._parent.delatt([att_name])


class FaceSet:
    """
    FaceSet
    """

    def __init__(self, filename: str, parent: "EltSet"):
        self.filename = filename
        self._parent = parent

    def __repr__(self):
        return str(self.filename)


class PSet:
    """
    Pset
    """

    def __init__(self, name: str, parent: PyLaGriT | MO):
        self.name = name
        self._parent = parent

    def __repr__(self):
        return str(self.name)

    def delete(self):
        cmd = "pset/" + self.name + "/delete"
        self._parent.sendcmd(cmd)

    def minmax_xyz(self):
        cmd = "/".join(
            [
                "cmo/printatt",
                cast(MO, self._parent).name,
                "-xyz-",
                "minmax",
                "pset,get," + self.name,
            ]
        )
        self._parent.sendcmd(cmd)

    def minmax(self, attname: Optional[str] = None, stride=(1, 0, 0)):
        cast(MO, self._parent).printatt(
            attname=attname, stride=stride, pset=self, ptype="minmax"
        )

    def list(self, attname=None, stride=(1, 0, 0)):
        cast(MO, self._parent).printatt(
            attname=attname, stride=stride, pset=self, ptype="list"
        )

    def setatt(self, attname: str, value: int | float | str):
        cmd = "/".join(
            [
                "cmo/setatt",
                cast(MO, self._parent).name,
                attname,
                "pset get " + self.name,
                str(value),
            ]
        )
        self._parent.sendcmd(cmd)

    def refine(
        self,
        refine_type="element",
        refine_option="constant",
        interpolation=" ",
        prange=(-1, 0, 0),
        field=" ",
        inclusive_flag="exclusive",
        prd_choice: Optional[int] = None,
    ):
        prange = [str(v) for v in prange]
        cmd = [
            "refine",
            refine_option,
            field,
            interpolation,
            refine_type,
            "pset get " + self.name,
            ",".join(prange),
            inclusive_flag,
        ]
        if prd_choice is not None:
            cmd.append("amr " + str(prd_choice))
        self._parent.sendcmd("/".join(cmd))

    def eltset(
        self,
        membership: Literal["inclusive", "exclusive", "face"] = "inclusive",
        name: Optional[str] = None,
    ):
        """
        Create eltset from pset

        Parameters
        ----------
        membership : str
            type of element membership, one of [inclusive,exclusive,face]

        name : str
            Name of element set to be used within LaGriT

        Returns
        -------
        EltSet

        """
        if name is None:
            name = new_name("e", cast(MO, self._parent).lg.core.eltset_names())
        cmd = ["eltset", name, membership, "pset", "get", self.name]
        self._parent.sendcmd("/".join(cmd))
        return EltSet(name, cast(MO, self._parent))

    def expand(
        self, membership: Literal["inclusive", "exclusive", "face"] = "inclusive"
    ):
        """
        Add points surrounding pset to pset

        membership : str
            type of element membership, one of [inclusive,exclusive,face]

        """
        e = self.eltset(membership=membership)
        self._parent.sendcmd("pset/" + self.name + "/delete")
        self = e.pset(name=self.name)

    def interpolate(
        self,
        method: Literal["map", "voronoi", "continuous", "default"],
        attsink: str,
        cmosrc: MO,
        attsrc: str,
        interp_function: Optional[
            Literal[
                "linear",
                "asinh",
                "log",
                "copy",
                "sequence",
                "min",
                "incmin",
                "max",
                "incmax",
                "and",
                "or",
            ]
            | int
            | float
        ] = None,
    ):
        """
        Interpolate values from attribute attsrc from mesh object cmosrc to current mesh object
        """
        cast(MO, self._parent).interpolate(
            method=method,
            attsink=attsink,
            stride=["pset", "get", self.name],
            cmosrc=cmosrc,
            attsrc=attsrc,
            interp_function=interp_function,
        )

    def interpolate_continuous(
        self,
        attsink,
        cmosrc,
        attsrc,
        interp_function=None,
        nearest: Optional[str] = None,
    ):
        cmd = [
            "intrp",
            "continuous",
            self.name + " " + attsink,
            ",".join(["pset", "get", self.name]),
            cmosrc.name + " " + attsrc,
        ]
        if nearest is not None:
            cmd += ["nearest", nearest]
        if interp_function is not None:
            cmd.append(interp_function)
        self._parent.sendcmd("/".join(cmd))

    def dump(self, filerootname: str, zonetype: Literal["zone", "zonn"] = "zone"):
        """
        Dump zone file of pset nodes

        Parameters
        ----------
        filerootname : str
            rootname of files to create, pset name will be added to name

        zonetype : str
            Type of zone file to dump, 'zone' or 'zonn'

        """
        cmd = ["pset", self.name, zonetype, filerootname + "_" + self.name, "ascii"]
        self._parent.sendcmd("/".join(cmd))

    def scale(
        self,
        scale_type: Literal["relative", "absolute"] = "relative",
        scale_geom: Literal["xyz", "rtz", "rtp"] = "xyz",
        scale_factor=(1, 1, 1),
        scale_center=(0, 0, 0),
    ):
        """
        Scale pset nodes by a relative or absolute amount

        Parameters
        ----------
        scale_type : str
            Scaling type may be 'relative' or 'absolute'

        scale_geom : str
            May be one of the geometry types 'xyz' (Cartesian), 'rtz' (cylindrical), or 'rtp' (spherical)

        scale_factor : Tuple [float, float, float]
            If scale_factor is relative, scaling factors are unitless multipliers. If absolute, scaling factors are constants added to existing coordinates.

        scale_center : Tuple [float, float, float]
            Geometric center to scale from

        """
        scale_factor = [str(v) for v in scale_factor]
        scale_center = [str(v) for v in scale_center]

        cmd = [
            "scale",
            ",".join(["pset", "get", self.name]),
            scale_type,
            scale_geom,
            ",".join(scale_factor),
            ",".join(scale_center),
        ]
        self._parent.sendcmd("/".join(cmd))

    def perturb(self, xfactor: float, yfactor: float, zfactor: float):
        """
        This command moves node coordinates in the following manner.

        Three pairs of random numbers between 0 and 1 are generated.
        These pairs refer to the x, y and z coordinates of the nodes respectively.
        The first random number of each pair is multiplied by the factor given in
        the command. The second random number is used to determine
        if the calculated offset is to be added or subtracted from the coordinate.
        """

        cmd = [
            "perturb",
            ",".join(["pset", "get", self.name]),
            str(xfactor),
            str(yfactor),
            str(zfactor),
        ]
        self._parent.sendcmd("/".join(cmd))

    def trans(self, xold: Tuple[float, float, float], xnew: Tuple[float, float, float]):
        """
        Translate points within a pset by the linear translation from (xold, yold, zold) to (xnew, ynew, znew)

        Parameters
        ----------
        xold : Tuple [float, float, float]
            Tuple containing point (xold, yold, zold) to translate from

        xnew : Tuple [float, float, float]
            Tuple containing point (xnew, ynew, znew) to

        """
        cmd = [
            "trans",
            ",".join(["pset", "get", self.name]),
            ",".join([str(v) for v in xold]),
            ",".join([str(v) for v in xnew]),
        ]
        self._parent.sendcmd("/".join(cmd))

    def smooth(self, *args, **kwargs):
        if "algorithm" not in kwargs:
            algorithm = " "
        else:
            algorithm = kwargs["algorithm"]
        cmd = ["smooth", "position", algorithm, "pset get " + self.name]
        for a in args:
            cmd.append(a)
        self._parent.sendcmd("/".join(cmd))

    def pset_attribute(
        self,
        attribute: str,
        value: int | float,
        comparison: Literal["lt", "le", "gt", "ge", "eq", "ne"] = "eq",
        name: Optional[str] = None,
    ):
        """
        Define PSet from another PSet by attribute

        Parameters
        ----------
        attribute : str
            Nodes defined by attribute name.

        value : int | float
            Attribute ID value.

        comparison : str
            Attribute comparison, default is eq.

        name : str
            The name to be assigned to the PSet created.

        Returns
        -------
        PSet

        """
        if name is None:
            name = new_name("p", cast(MO, self._parent).lg.core.pset_names())

        cmd = "/".join(
            [
                "pset",
                name,
                "attribute " + attribute,
                "pset,get," + self.name,
                " " + comparison + " " + str(value),
            ]
        )

        self._parent.sendcmd(cmd)
        return PSet(name, self._parent)
