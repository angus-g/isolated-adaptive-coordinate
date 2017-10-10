from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

npy_include_dir = numpy.get_include()

ext_modules = [Extension('adaptive', ['adaptive/adaptive_wrap.pyx'],
                         libraries=['gfortran'],
                         include_dirs=[npy_include_dir],
                         extra_objects=['adaptive/adaptive.o', 'adaptive/MOM_EOS.o', 'adaptive/MOM_EOS_Wright.o', 'adaptive/MOM_EOS_linear.o']),
               Extension('remapping', ['remapping/remapping_wrap.pyx'],
                         libraries=['gfortran'],
                         include_dirs=[npy_include_dir],
                         extra_objects=['remapping/remapping.o', 'remapping/MOM_remapping.o', 'remapping/MOM_error_handler.o',
                                        'remapping/PCM_functions.o', 'remapping/PLM_functions.o', 'remapping/PPM_functions.o',
                                        'remapping/PQM_functions.o', 'remapping/polynomial_functions.o', 'remapping/regrid_solvers.o',
                                        'remapping/regrid_edge_slopes.o', 'remapping/regrid_edge_values.o'])]

setup(name='Adaptive Coordinate',
      cmdclass={'build_ext': build_ext},
      ext_modules=ext_modules)
