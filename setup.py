# Cython compile instructions
import numpy
from setuptools import setup, Extension
from Cython.Build import build_ext

# Use python setup.py build --inplace
# to compile

extensions = [
    Extension("pyshot",
              sources=["pyshot.pyx", './src/shot_descriptor.cpp'],
              include_dirs=[
                  numpy.get_include(),
                  'include/',
                  '/usr/include/eigen3/'],
              libraries=["lz4"],
              extra_compile_args=["-O3"],
              extra_link_args=['-L/usr/include/'],
              language="c++",
              define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
)
]

setup(
  name="pyshot",
  cmdclass={'build_ext': build_ext},
  ext_modules=extensions,
)