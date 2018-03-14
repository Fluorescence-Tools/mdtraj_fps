from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

extension = Extension(
        '_fps',
        sources=['./mdtraj_fps/_fps.pyx', './md_traj/mt19937cok.cpp'],
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
    cmdclass={'build_ext': build_ext},
    ext_modules=[extension]
)