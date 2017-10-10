adaptive.cpython-36m-x86_64-linux-gnu.so: adaptive_wrap.pyx setup.py adaptive.o
	python setup.py build_ext -i

# Fortran dependencies
adaptive.o: MOM_EOS.o
MOM_EOS.o: MOM_EOS_linear.o MOM_EOS_Wright.o

%.o: %.F90
	gfortran -fcray-pointer -fdefault-double-8 -fdefault-real-8 -Waliasing -ffree-line-length-none -g -Wall -fPIC -c $<

clean:
	rm -f *.o *.mod *.so adaptive_wrap.c
