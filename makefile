all: main

COMPILER = mpifort
COMPFLAG = -Wall -Ofast
PETSC = -I${PETSC_DIR}/include -I${PETSC_DIR}/arch-linux2-c-debug/include

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

objects = props.o ptcl_props.o petsc_lib.o elec_lib.o ion_lib.o lapl_lib.o eqn_lib.o

#--------------------------------------------------------------------------
debug: COMPFLAG += -Wall -Wextra -pedantic -g -O0 -fimplicit-none -fbacktrace -fcheck=all
debug: main
#--------------------------------------------------------------------------

main: $(objects) main.o chkopts
	-${FLINKER} $(COMPFLAG) -o main $(objects) main.o ${PETSC_KSP_LIB}

%.o: %.F90
	$(COMPILER) $(COMPFLAG) -c $< $(PETSC) 
