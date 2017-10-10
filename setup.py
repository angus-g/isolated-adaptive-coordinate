from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

npy_include_dir = numpy.get_include()

ext_modules = [Extension("adaptive", ["adaptive_wrap.pyx"],
                         libraries=['gfortran'],
                         include_dirs = [npy_include_dir],
                         extra_objects=["adaptive.o", "MOM_EOS.o", "MOM_EOS_Wright.o", "MOM_EOS_linear.o"])]

setup(name='Adaptive Coordinate',
      cmdclass={'build_ext': build_ext},
      ext_modules=ext_modules)
