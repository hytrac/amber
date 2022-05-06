# Compiler
FC     = ifort
FF     = -align -fpp -ftz -g -ip -ipo -O3 -traceback
DEBUG  = -check -debug
INTEL  = -qmkl=parallel -mcmodel=medium -shared-intel
OPENMP = -qopenmp
FLAGS  = $(FF) $(INTEL) $(OPENMP)
#FLAGS  = $(FF) $(DEBUG) $(INTEL) $(OPENMP)

# INC and LIB
INC = 
LIB = 

# Object files
OBJ = constant.o timing.o mkl.o helper.o input.o \
      sim.o cosmo.o mesh.o esf.o grf.o lpt.o reion.o cmb.o h21cm.o \
      simulation.o cosmology.o meshmake.o esfhalo.o grfmake.o lptmake.o \
      reionization.o cmbreion.o h21cmreion.o

# Executable
amber.x : $(OBJ)
	$(FC) $(FLAGS) $(OBJ) $(INC) $(LIB) amber.f90 -o $@
	rm *.o *.mod

%.o: %.f90
	$(FC) $(FLAGS) $(INC) -c $*.f90

all : amber.x

clean :
	rm *.o *.mod *.x
