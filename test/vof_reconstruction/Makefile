FC = mpif90
FOPTS = -g -Wall -fcheck=bound
CODE_DIR = $(REO)

#HYPRE_DIR =
INCLUDES = -I$(HYPRE_DIR)/include -I$(CODE_DIR)
LIBS = -L$(HYPRE_DIR)/lib -lHYPRE -lm
FFLAGS = $(FOPTS) $(INCLUDES)
FLAGS = $(OPTS) $(LIBS) -lstdc++

#===============================================================================
# TEST FOR POISSON EQUATION
#===============================================================================
reconstruction: reconstruction.o
	$(FC) -o $@ $< grid_2d.o volume_of_fluid.o boundary_conditions.o $(FLAGS)

reconstruction.o: reconstruction.f08 grid_2d.o volume_of_fluid.o boundary_conditions.o
	$(FC) $(FFLAGS) -c $<
#===============================================================================

#===============================================================================
# Code
#===============================================================================
volume_of_fluid.o: $(CODE_DIR)/volume_of_fluid.f08 boundary_conditions.o grid_2d.o
	$(FC) $(FFLAGS) -c $<

boundary_conditions.o: $(CODE_DIR)/boundary_conditions.f08 grid_2d.o
	$(FC) $(FFLAGS) -c $<

grid_2d.o: $(CODE_DIR)/grid_2d.f08
	$(FC) $(FFLAGS) -c $<

plot:
	gnuplot plot.pl

clean:
	rm -r *.o *.mod reconstruction h plot.png
