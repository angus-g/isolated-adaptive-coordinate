FFLAGS ?= -fcray-pointer -fdefault-double-8 -fdefault-real-8 -Waliasing -ffree-line-length-none -g -Wall -fPIC

# Wrapped adaptive module
adaptive.cpython-36m-x86_64-linux-gnu.so remapping.cpython-36m-x86_64-linux-gnu.so: adaptive/adaptive_wrap.pyx setup.py adaptive/adaptive.o remapping/remapping_wrap.pyx remapping/remapping.o
	python setup.py build_ext -i

# Fortran dependencies for adaptive module
adaptive/adaptive.o: adaptive/MOM_EOS.o
adaptive/MOM_EOS.o: adaptive/MOM_EOS_linear.o adaptive/MOM_EOS_Wright.o

# Fortran dependencies for remapping module
remapping/remapping.o: remapping/MOM_remapping.o
remapping/MOM_remapping.o: remapping/MOM_error_handler.o remapping/regrid_edge_values.o remapping/regrid_edge_slopes.o remapping/PCM_functions.o remapping/PLM_functions.o remapping/PPM_functions.o remapping/PQM_functions.o
remapping/PPM_functions.o: remapping/regrid_edge_values.o
remapping/PQM_functions.o: remapping/regrid_edge_values.o
remapping/regrid_edge_values.o: remapping/regrid_solvers.o remapping/polynomial_functions.o
remapping/regrid_edge_slopes.o: remapping/regrid_solvers.o remapping/polynomial_functions.o
remapping/regrid_solvers.o: remapping/MOM_error_handler.o

%.o: %.F90
	gfortran $(FFLAGS) -c -J $(dir $@) -o $@ $<

clean:
	rm -f adaptive/*.o adaptive/*.mod adaptive/adaptive_wrap.c remapping/*.o remapping/*.mod remapping/remapping_wrap.c *.so
