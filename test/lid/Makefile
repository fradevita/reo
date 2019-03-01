FC = mpif90
FOPTS = -g -Wall -O2 -fcheck=bound -fdefault-real-8
CODE_DIR = $(REO)

#HYPRE_DIR = 
INCLUDES = -I$(HYPRE_DIR)/include -I$(CODE_DIR)
LIBS = -L$(HYPRE_DIR)/lib -lHYPRE -lm
FFLAGS = $(FOPTS) $(INCLUDES)
FLAGS = $(OPTS) $(LIBS) -lstdc++

#===============================================================================
# TEST FOR POISSON EQUATION
#===============================================================================
lid: lid.o
	$(FC) -o $@ $< grid_2d.o boundary_conditions.o navier_stokes.o hypre_solver.o $(FLAGS)

lid.o: lid.f08 grid_2d.o navier_stokes.o
	$(FC) $(FFLAGS) -c $<
#===============================================================================

#===============================================================================
# Code
#===============================================================================
navier_stokes.o: $(CODE_DIR)/navier_stokes.f08 grid_2d.o boundary_conditions.o hypre_solver.o
	$(FC) $(FFLAGS) -c $<

hypre_solver.o: $(CODE_DIR)/hypre_solver.f08 grid_2d.o boundary_conditions.o
	$(FC) $(FFLAGS) -c $<

boundary_conditions.o: $(CODE_DIR)/boundary_conditions.f08 grid_2d.o
	$(FC) $(FFLAGS) -c $<

grid_2d.o: $(CODE_DIR)/grid_2d.f08
	$(FC) $(FFLAGS) -c $<

plot:
	gnuplot plot.pl

clean:
	rm -r *.o *.mod lid log *.png out xprof yprof
