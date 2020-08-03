from setuptools import setup
from Cython.Build import cythonize

setup(ext_modules=cythonize("1D_TISE.pyx",annotate=True))
