from . import _fps
import ctypes as C
import platform
import os
import json
import numpy as np


b, o = platform.architecture()
package_directory = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(package_directory, './dll')

if 'Windows' in o:
    if '32' in b:
        fpslibrary = 'fpsnative.win32.dll'
    elif '64' in b:
        fpslibrary = 'fpsnative.win64.dll'
else:
    if platform.system() == 'Linux':
        fpslibrary = 'liblinux_fps.linux32.so'
    else:
        fpslibrary = 'libav.dylib'

fps_dll_file = os.path.join(path, fpslibrary)
_fps_dll = np.ctypeslib.load_library(fps_dll_file, './dll')

_fps_dll.calculate1R.restype = C.c_int
_fps_dll.calculate1R.argtypes = [C.c_double, C.c_double, C.c_double,
                                 C.c_int, C.c_double,
                                 C.POINTER(C.c_double), C.POINTER(C.c_double), C.POINTER(C.c_double),
                                 C.POINTER(C.c_double), C.c_int, C.c_double,
                                 C.c_double, C.c_int,
                                 C.POINTER(C.c_char)]

_fps_dll.calculate3R.restype = C.c_int
_fps_dll.calculate3R.argtypes = [C.c_double, C.c_double, C.c_double, C.c_double, C.c_double,
                                 C.c_int, C.c_double,
                                 C.POINTER(C.c_double), C.POINTER(C.c_double), C.POINTER(C.c_double),
                                 C.POINTER(C.c_double), C.c_int, C.c_double,
                                 C.c_double, C.c_int,
                                 C.POINTER(C.c_char)]

_periodic_table = json.load(open(os.path.join(package_directory, 'elements.json')))['Periodic Table']
VDW_DICT = dict((key, _periodic_table[key]["vdW radius"])
                for key in _periodic_table.keys())


def get_vdw(trajectory):
    """Get a vector of the vdw-radii
    :param trajectory: mdtraj
    :return:
    """
    return np.array([_periodic_table[atom.element.symbol]['vdW radius'] for atom in trajectory.topology.atoms],
                    dtype=np.float64)


def pRDA(distances=None, av1=None, av2=None, bins=None, n_samples=10000, normed=True):
    """Calculates the distance distribution between two accessible volumes
    """
    if isinstance(distances, np.ndarray):
        d = distances
    else:
        d = random_distances(av1, av2, n_samples)
    if bins is None:
        bins = np.linspace(10, 100, 90)
    return np.histogram(d, bins=bins, normed=normed)


def RDAMean(distances=None, av1=None, av2=None, nSamples=50000):
    """Calculate the mean distance between two accessible volumes
    """
    if isinstance(distances, np.ndarray):
        return distances.mean()
    else:
        return _fps.RDAMean(av1, av2, nSamples)


def RDAMeanE(forster_radius, distances=None, av1=None, av2=None, nSamples=50000):
    """Calculate the FRET-averaged (PDA/Intensity) distance between two accessible volumes
    """
    if isinstance(distances, np.ndarray):
        transfer = (1./(1.+(distances/forster_radius)**6.0))
        mean_transfer = transfer.mean()
        return (1./mean_transfer - 1.)**(1./6.) * forster_radius
    else:
        return _fps.RDAMeanE(av1, av2, forster_radius, nSamples)


def dRmp(av1, av2):
    """Calculate the distance between the mean position of two accessible volumes
    """
    return np.sqrt(((av1.Rmp - av2.Rmp) ** 2).sum())


def random_distances(av1, av2, n_samples=10000):
    """Generates a set of random distances
    """
    return _fps.random_distances(av1, av2, n_samples)


def density2points(n, npm, dg, density, r0, ng):
    """
    :param n:
    :param npm:
    :param dg:
    :param density:
    :param r0:
    :param ng:
    :return:
    """
    return _fps.density2points(n, npm, dg, density, r0, ng)


def calculate_1_radius(x, y, z, vdw, l, w, r1, atom_i, linkersphere=0.5, linknodes=3, vdwRMax=1.8, dg=0.5, **kwargs):
    """
    :param l: float
        linker length
    :param w: float
        linker width
    :param r: float
        dye-radius
    :param atom_i: int
        attachment-atom index
    :param x: array
        Cartesian coordinates of atoms (x) in angstrom
    :param y: array
        Cartesian coordinates of atoms (y) in angstrom
    :param z: array
        Cartesian coordinates of atoms (z) in angstrom
    :param vdw:
        Van der Waals radii (same length as number of atoms)
    :param linkersphere: float
        Initial linker-sphere to start search of allowed dye positions
    :param linknodes: int
        By default 3
    :param vdwRMax: float
        Maximal Van der Waals radius
    :param dg: float
        Resolution of accessible volume in Angstrom
    :param verbose: bool
        If true informative output is printed on std-out


    """
    n_atoms = len(vdw)

    npm = int(np.floor(l / dg))
    ng = 2 * npm + 1
    ng3 = ng * ng * ng
    density = np.zeros(ng3, dtype=np.uint8)
    x0, y0, z0 = x[atom_i], y[atom_i], z[atom_i]
    r0 = np.array([x0, y0, z0], dtype=np.float64)
    vdw = vdw.astype(np.float64, order='C')

    _x = x.ctypes.data_as(C.POINTER(C.c_double))
    _y = y.ctypes.data_as(C.POINTER(C.c_double))
    _z = z.ctypes.data_as(C.POINTER(C.c_double))
    _vdw = vdw.ctypes.data_as(C.POINTER(C.c_double))

    _density = density.ctypes.data_as(C.POINTER(C.c_char))
    n = _fps_dll.calculate1R(l, w, r1, atom_i, dg, _x, _y, _z, _vdw, n_atoms, vdwRMax, linkersphere, linknodes,
                             _density)

    density = density.astype(np.float32)
    points = density2points(n, npm, dg, density, r0, ng)
    density = density.reshape([ng, ng, ng])
    return points, density, ng, r0


def calculate_3_radius(x, y, z, vdw, l, w, r1, r2, r3, atom_i, linkersphere=0.5, linknodes=3, vdwRMax=1.8, dg=0.4,
                       **kwargs):
    """
    :param l: float
        linker length
    :param w: float
        linker width
    :param r1: float
        Dye-radius 1
    :param r2: float
        Dye-radius 2
    :param r3: float
        Dye-radius 3
    :param atom_i: int
        attachment-atom index
    :param x: array
        Cartesian coordinates of atoms (x)
    :param y: array
        Cartesian coordinates of atoms (y)
    :param z: array
        Cartesian coordinates of atoms (z)
    :param vdw:
        Van der Waals radii (same length as number of atoms)
    :param linkersphere: float
        Initial linker-sphere to start search of allowed dye positions
    :param linknodes: int
        By default 3
    :param vdwRMax: float
        Maximal Van der Waals radius
    :param dg: float
        Resolution of accessible volume in Angstrom
    :param verbose: bool
        If true informative output is printed on std-out
    :return:

    """
    n_atoms = len(vdw)

    npm = int(np.floor(l / dg))
    ng = 2 * npm + 1
    ng3 = ng * ng * ng
    density = np.zeros(ng3, dtype=np.uint8)
    x0, y0, z0 = x[atom_i], y[atom_i], z[atom_i]
    r0 = np.array([x0, y0, z0], dtype=np.float64)

    _x = x.ctypes.data_as(C.POINTER(C.c_double))
    _y = y.ctypes.data_as(C.POINTER(C.c_double))
    _z = z.ctypes.data_as(C.POINTER(C.c_double))
    _vdw = vdw.ctypes.data_as(C.POINTER(C.c_double))

    _density = density.ctypes.data_as(C.POINTER(C.c_char))
    n = _fps_dll.calculate3R(l, w, r1, r2, r3, atom_i, dg, _x, _y, _z, _vdw, n_atoms, vdwRMax, linkersphere, linknodes,
                             _density)

    density = density.astype(np.float32)
    points = density2points(n, npm, dg, density, r0, ng)
    density = density.reshape([ng, ng, ng])
    return points, density, ng, r0