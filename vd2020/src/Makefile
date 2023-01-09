#PROJROOT := $(shell pwd | sed 's,source,,g')

#BUILDDIR = $(PROJROOT)/build
#BINDIR = $(PROJROOT)/bin

SRCDIR = .
LLTB = ../wLLTBBackground
BINDIR = ../wLLTBBackground/bin
BUILDDIR = ../wLLTBBackground/build
DEPSDIR = $(BUILDDIR)/deps


F90C ?= gfortran
F90FLAGS ?= -O3 -ffree-form -ffree-line-length-none -fopenmp -std=gnu -fbacktrace -g -cpp -fPIC

INCLUDE = -I$(BUILDDIR) -I$(LLTB)/include
LIBS = -L$(BUILDDIR) -L$(LLTB)/lib -lwlltbpack

ifeq ($(F90C),mpif90)
    thiscompiler := $(shell  mpif90 --showme -v | awk '{print $1}')
else
    thiscompiler := $(F90C)
endif
ifeq ($(thiscompiler),gfortran)
MODULEFLAG = -J
endif
ifeq ($(thiscompiler),ifort)
MODULEFLAG = -module
endif
ifeq ($(MODULEFLAG),)
MODULEFLAG = -J
endif
MODULES = $(MODULEFLAG) $(BUILDDIR)

PLAINFOBJ := $(shell ls *f90)
FOBJ = $(PLAINFOBJ:%.f90=$(BUILDDIR)/%.o)
OBJ = $(FOBJ)




default: driver


$(BUILDDIR)/vdii_tools.o: $(BUILDDIR)/vdii_constants.o
$(BUILDDIR)/voidnrodeint.o: $(BUILDDIR)/voidnrtype.o
$(BUILDDIR)/voiddistances.o: $(BUILDDIR)/void_dvodeinterface.o $(BUILDDIR)/vdii_vdipatch_FakeCMBParams.o $(BUILDDIR)/voidnrodeint.o $(BUILDDIR)/precision.o $(BUILDDIR)/vdii_tools.o
$(BUILDDIR)/voiddistances.o:
$(BUILDDIR)/driver.o: $(BUILDDIR)/voiddistances.o

driver: $(OBJ)
	$(F90C) -o $(BINDIR)/vd2020 $(OBJ) $(LIBS) $(INCLUDE) $(F90FLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.f90
	$(F90C) $(F90FLAGS) $(MODULES) $(INCLUDE)  -c $< -o $@

# Special treatment: external library, which works, so suppress warnings.
$(BUILDDIR)/void_dvodeinterface.o: $(SRCDIR)/void_dvodeinterface.f90
	$(F90C) -w $(F90FLAGS) $(MODULES) $(INCLUDE) -c $< -o $@

clean:
	rm -f $(BUILDDIR)/*.mod $(BUILDDIR)/*.o
