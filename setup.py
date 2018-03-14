from setuptools import setup
from setuptools.extension import Extension
from Cython.Distutils import build_ext
import numpy

extension = Extension(
        'mdtraj_fps',
        sources=['./mdtraj_fps/_fps.pyx', './mdtraj_fps/mt19937cok.cpp'],
        include_dirs=[numpy.get_include(), "."],
        # extra_compile_args=['-fopenmp'],
        #extra_link_args=['-fopenmp'],
        libraries=[],
        language="c++"
)

extension.pyrex_directives = {"boundscheck": False,
                              "wraparound": False,
                              "cdivision": True,
                              "profile": False
}

setup(
    name="mdtraj_fps",
    version="0.0.1",
    author="Thomas Peulen",
    author_email="thomas.otavio.peulen@gmail.com",
    description=("A tool to calculate FRET observables for MD trajectories by accessible volume calculations."
                 "Here, spatial density of flexible coupled dyes/labels is approximated by the sterically allowed"
                 "space for labels modeled by ellipsoids attached by a flexible cylinder."
                 ),
    license="LGPL",
    cmdclass={'build_ext': build_ext},
    ext_modules=[extension]
)