SHELL = /bin/bash

PROJROOT := $(shell pwd | sed 's,source,,g')

SRCDIR = $(PROJROOT)/source
BUILDDIR = $(PROJROOT)/build
DEPSDIR = $(PROJROOT)/build/deps
BINDIR = $(PROJROOT)/bin
INCDIR = $(PROJROOT)/include
LIBDIR = $(PROJROOT)/lib

CHKBUILD := $(shell mkdir -p $(BUILDDIR); mkdir -p $(DEPSDIR); mkdir -p $(BINDIR); mkdir -p $(INCDIR); mkdir -p $(LIBDIR))

LLTB = $(PROJROOT)/../LLTBBackground

F90C = gfortran
F90FLAGS = -fbacktrace -g -O3 -ffree-form -shared -ffree-line-length-none -fopenmp -cpp -fPIC
#F90FLAGS = -fbacktrace -g -O3 -ffree-form -ffree-line-length-none -fopenmp -Wl,-framework -Wl,CoreFoundation -cpp
#F90FLAGS = -fbacktrace -g -O3 -ffree-form -ffree-line-length-none -fopenmp  -cpp
INCLUDE = -I$(BUILDDIR)
MODULEFLAG = -J
MODULES = $(MODULEFLAG) $(BUILDDIR)
LIBS = -I$(LLTB)/include -L$(LLTB)/lib -llltbpack

#ifneq ($(HOME),/Users/Valkenburg)
have_a := $(shell if [[ -f Makefile.local ]]; then echo LOCAL; else echo NONLOCAL; fi) 
#else
#have_a := NONLOCAL
#endif

ifeq ($(strip $(have_a)),LOCAL) 
include Makefile.local
endif


OBJ = \
	constants.o \
	wlltb_rombint.o \
	kofrattributes.o \
	types.o \
	wlltb_types_funcs.o \
	wlltb_errors.o \
	kofr.o \
	wlltb_funcs.o \
	wlltb_dvodeinterface.o \
	wlltb_integratormod.o \
	wlltb_background.o \
	driver.o

INCOBJ = $(OBJ:%=$(BUILDDIR)/%)
#F90FLAGS +=  $(OLDCAMB) $(MODULES)
LLFFLAGS = $(F90FLAGS) $(OLDCAMB) $(MODULES)

export FFLAGS
export F90C

GETDEPENDENCE = $(strip $(foreach a, $2, $(if $(findstring $(a), $1),  $(MYVAR), $(eval MYVAR = $a))))
MAKEDEPENDENCE = mkdir -p $(DEPSDIR); echo $@: $(call GETDEPENDENCE, $@, $(INCOBJ)) > $(@:$(BUILDDIR)/%.o=$(DEPSDIR)/%.d)

ifeq ($(strip $(have_a)),LOCAL) 

default: message binclean driver

message:
	@echo ""; echo ""; echo "         Makefile.local found. Using that in stead of Makefile."; echo ""; echo "";
	@echo $(HOME)
else

default: binclean driver wlltbpack

endif 

all: clean driver wlltbpack

-include $(OBJ:%.o=$(DEPSDIR)/%.d)

$(BUILDDIR)/wlltb_dvodeinterface.o: $(SRCDIR)/dvode/dvode_f90_m.f90

$(BUILDDIR)/%.o: $(SRCDIR)/%.f90
	$(F90C) $(LLFFLAGS) $(INCLUDE) -c $< -o $@ $(LIBS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.F90
	$(F90C) $(LLFFLAGS) $(INCLUDE) -c $< -o $@ $(LIBS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CINCLUDE) -c $*.c

driver: $(INCOBJ)
	$(F90C) $(LLFFLAGS) $() -o $(BINDIR)/pert  $(INCOBJ) $(LIBS)
	
clean:
	rm -f $(BUILDDIR)/*o $(BUILDDIR)/*mod $(BINDIR)/pert $(LIBDIR)/*.a $(INCDIR)/*.mod
	
binclean:
	rm -f $(BINDIR)/pert


# TESTING DVODE:
$(BUILDDIR)/dvode_f90_m.o: dvode/dvode_f90_m.f90
	$(F90C) $(LLFFLAGS) $(INCLUDE) -c dvode/dvode_f90_m.f90 -o $(BUILDDIR)/dvode_f90_m.o 

dvode: $(BUILDDIR)/dvode_f90_m.o $(BINDIR)/dvode

$(BINDIR)/dvode:
	$(F90C)  -std=legacy  $(LLFFLAGS) $(INCLUDE) -o $(BINDIR)/dvode dvode/example1.f90 $(BUILDDIR)/dvode_f90_m.o 
	
wlltbpack: $(LLTB)/lib/liblltbpack.a $(INCOBJ)
	mkdir -p $(BUILDDIR)/tmppackage; touch $(BUILDDIR)/tmppackage/abc;	rm -rf $(BUILDDIR)/tmppackage/*;	#cp $(LIBDIR)/*a $(BUILDDIR)/tmppackage/.
	cp $(LLTB)/lib/liblltbpack.a $(BUILDDIR)/tmppackage/.
	cd $(BUILDDIR)/tmppackage &&  rm -f libwlltbpack.a &&	ar -x liblltbpack.a &&	ar -csrv libwlltbpack.a ../*.o *.o
	mv $(BUILDDIR)/tmppackage/libwlltbpack.a $(LIBDIR)/.
	rm -rf $(BUILDDIR)/tmppackage
	cp $(BUILDDIR)/*.mod $(INCDIR)/.
	cp $(LLTB)/include/*.mod $(INCDIR)/.
	
$(LLTB)/lib/liblltbpack.a:
	cd $(LLTB)/src && make lltbpackage F90C="$(F90C)" FFLAGS="$(F90FLAGS)"
	
#  	
