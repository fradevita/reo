FC = mpif90
FOPTS = -g -Wall -O2 -fcheck=bound -fdefault-real-8

#HYPRE_DIR = /scratch/fdv/hypre-2.11.2/src/hypre
INCLUDES = -I$(HYPRE_DIR)/include
LIBS = -L$(HYPRE_DIR)/lib -lHYPRE -lm
FFLAGS = $(FOPTS) $(INCLUDES)
FLAGS = $(OPTS) $(LIBS) -lstdc++

#===============================================================================
# TESTING DURING CODING
#===============================================================================
testing: testing.o
	$(FC) -o $@ $< grid_2d.o navier_stokes.o boundary_conditions.o hypre_solver.o $(FLAGS)

testing.o: testing.f08 navier_stokes.o grid_2d.o
	$(FC) $(FFLAGS) -c $<
#===============================================================================

#===============================================================================
# Code
#===============================================================================
navier_stokes.o: navier_stokes.f08 grid_2d.o boundary_conditions.o hypre_solver.o
	$(FC) $(FFLAGS) -c $<

hypre_solver.o: hypre_solver.f08 grid_2d.o boundary_conditions.o
	$(FC) $(FFLAGS) -c $<

boundary_conditions.o: boundary_conditions.f08 grid_2d.o
	$(FC) $(FFLAGS) -c $<

grid_2d.o: grid_2d.f08
	$(FC) $(FFLAGS) -c $<

clean:
	rm -r *.o *.mod testing fort.*
