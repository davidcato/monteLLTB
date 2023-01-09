# ============================================================================
# Name        : EllipticLTB
# Author      : 
# Version     :
# Copyright   : Your copyright notice
# Description : Makefile for Hello World in Fortran
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
#FORTRAN_COMPILER = g95
# EllipticLTB with test driver.
ELLIPSI = ../../EllipticIntegrals

# Fortran Multiple Precision package by D.M. Smith:
# http://myweb.lmu.edu/dmsmith/FMLIB.html
FMLIB = 
#../../FMLib

IFLAG = -I

#ifneq ($(HOME),/Users/Valkenburg)
#F90C     = ifort
#FFLAGS = -openmp -O2 -ip -W0 -WB -fpp2 -vec_report0 -i-static -traceback -g
F90C = gfortran
FFLAGS = -fbacktrace -g -O3 -ffree-form -ffree-line-length-none -fopenmp -Wl,-framework -Wl,CoreFoundation -cpp -fPIC
DEBUGFLAGS = -check all -check noarg_temp_created -debug-parameters all
#-fp-stack- -fp-model precise 
#F90C     = ifort
#FFLAGS = -openmp -O2 -ip -W0 -WB -fpp2 -vec_report0 -i-static -traceback -g 
#else
#F90C = gfortran
#FFLAGS = -fbacktrace -g -O3 -ffree-form -ffree-line-length-none -fopenmp -Wl,-framework -Wl,CoreFoundation -cpp -std=gnu 
#DEBUGFLAGS = -check all -check noarg_temp_created -debug-parameters all -fp-stack-check -warn all -implicitnone  -fp-model precise
#endif

#FFLAGS += ${DEBUGFLAGS}

#F90C   = g95
#FFLAGS = -ftrace=full -O3 -ffree-form -ffree-line-length-huge -w -Wprecision-loss

INCLUDES = -I$(ELLIPSI)/include

DRIVER        = driver.f90

LIBS =  -L$(ELLIPSI)/lib -lelliptics 
LIBS += -Wl,-framework,CoreFoundation,-stack_size,0x10000000 -I/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.10.sdk/usr/include/ -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.10.sdk/usr/lib/

ifneq ($(FMLIB),)
INCLUDES += -I$(FMLIB)/include
LIBS += -L$(FMLIB)/lib -lfm
FFLAGS += -DFMLIB
endif

F90FLAGSLOCAL      = $(FFLAGS) $(INCLUDES)

FC            = $(F90C)

OBJ       = lltb_precision.o lltb_complex_tools.o lltb_apprec.o lltb_params.o \
		lltb_background_errors.o \
		lltb_background_eqs.o lltb_prefuncs.o lltb_background_ht_ddok_funcs.o \
		lltb_background_ht_ddok.o lltb_background_ht_funcs.o lltb_background_ht.o \
		lltb_background_aofht_funcs.o lltb_background_aofht.o \
		lltb_background_limits.o\
		lltb_background.o

# default: clean test lltb lltbpackage
default: test lltb lltbpackage

all: $(ELLIPSI)/lib/libelliptics.a lltb

lltb_background.o: lltb_ricci.f90

%.o: %.f90
	$(F90C) $(F90FLAGSLOCAL) -c $*.f90 

%.o: %.F90
	$(F90C) $(F90FLAGSLOCAL) -c $*.F90 

test: $(OBJ) $(DRIVER) 
	mkdir -p ../bin
	$(F90C) $(DRIVER) -o ../bin/$@ $(OBJ) $(LIBS) $(F90FLAGSLOCAL)

lltb: $(ELLIPSI)/lib/libelliptics.a $(OBJ)
	mkdir -p ../lib
	mkdir -p ../include
	ar -csrv liblltb.a ${OBJ} 
	mv liblltb.a ../lib/.
	cp *.mod ../include/.
	cp $(ELLIPSI)/include/* ../include/.
	cp $(ELLIPSI)/lib/* ../lib/.
ifneq ($(FMLIB),)
	cp $(FMLIB)/lib/* ../lib/.
	cp $(FMLIB)/include/* ../include/.
endif

lltbpackage: lltb
	./lltbpackage.sh

$(ELLIPSI)/lib/libelliptics.a:
	cd ${ELLIPSI}/src && make all F90C="$(F90C)" FFLAGS="$(FFLAGS)"

clean:
	-rm -f *.o *.d *.a core *.mod test ../bin/test ../include/* ../lib/*
  
allclean: clean
	cd $(ELLIPSI)/src && make clean
