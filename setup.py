from setuptools import setup
from Cython.Build import cythonize

setup(ext_modules=cythonize("1D_TISE_c.pyx",annotate=True))
