# ============================================================================
# Name        : EllipticLTB
# Author      : 
# Version     :
# Copyright   : Your copyright notice
# Description : Makefile for Hello World in Fortran
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
#FORTRAN_COMPILER = gfortran
# EllipticLTB with test driver.

IFLAG = -I

#ifneq ($(HOME),/Users/Valkenburg)
#F90C     = ifort
#FFLAGS = -openmp -O2 -ip -W0 -WB -fpp2 -vec_report0 -i-static -traceback -g
F90C = gfortran
#FFLAGS = -fbacktrace -shared -g -O3 -ffree-form -ffree-line-length-none -fopenmp -Wl,-framework -Wl,CoreFoundation -cpp -fPIC
FFLAGS = -fbacktrace -shared -g -O3 -ffree-form -ffree-line-length-none -fopenmp -cpp -fPIC
DEBUGFLAGS = -check all -check noarg_temp_created -fp-model precise -debug-parameters all -fp-stack-
#F90C     = ifort
#FFLAGS = -openmp -O2 -ip -W0 -WB -fpp2 -vec_report0 -i-static -traceback -g 
#else
#F90C = gfortran
#FFLAGS = -fbacktrace -g -O3 -ffree-form -ffree-line-length-none -fopenmp -Wl,-framework -Wl,CoreFoundation -cpp -std=gnu 
#DEBUGFLAGS = -check all -check noarg_temp_created -fp-model precise -debug-parameters all -fp-stack-
#endif

#F90C   = g95
#FFLAGS = -ftrace=full -O3 -ffree-form -ffree-line-length-huge 

#INCLUDES = -I$(HEALPIX)/include -I$(cfitsio)/include -I$(MPI_OPENMPI_PATH)/include/
INCLUDES =

DRIVER        = driver.f90

F90FLAGS      = $(FFLAGS) $(INCLUDES)

LIBS = 
#-L$(HEALPIX)/lib -lhealpix -L$(cfitsio)/lib -lcfitsio -L$(MPI_OPENMPI_PATH)/lib -lmpi
LIBS += -Wl,-framework #,CoreFoundation ,-stack_size,0x10000000 -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.10.sdk/usr/include/ -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.10.sdk/usr/lib/

FC            = $(F90C)

OBJ       =  elliptics.o
#mympi.o cheese_types.o swiss_precision.o utils.o ch_functions.o geo_functions.o NR_odeint_geo.o NR_odeint_hl.o hl_functions.o swiss_geodesics.o inifile.o cheese_misc.o simpledistance.o exactdistance.o

default: clean test ellipticlib

all: ellipticlib

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90

test: $(OBJ) $(DRIVER)
	mkdir -p ../bin
	$(F90C) $(F90FLAGS) $(OBJ) $(DRIVER) -o ../bin/$@ $(LIBS)

ellipticlib: $(OBJ)
	mkdir -p ../lib
	mkdir -p ../include
	ar -csrv libelliptics.a elliptics.o 
	mv libelliptics.a ../lib/.
	cp *.mod ../include/.

clean:
	-rm -f *.o *.a *.d core *.mod test ../bin/* ../lib/* ../include/*
