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
poisson_dirichlet: poisson_dirichlet.o
	$(FC) -o $@ $< grid_2d.o hypre_solver.o $(FLAGS)

poisson_dirichlet.o: poisson_dirichlet.f08 grid_2d.o hypre_solver.o
	$(FC) $(FFLAGS) -c $<
#===============================================================================

#===============================================================================
# Code
#===============================================================================
hypre_solver.o: $(CODE_DIR)/hypre_solver.f08 grid_2d.o
	$(FC) $(FFLAGS) -c $<

grid_2d.o: $(CODE_DIR)/grid_2d.f08
	$(FC) $(FFLAGS) -c $<

plot:
	gnuplot plot.pl

clean:
	rm -r *.o *.mod poisson_dirichlet error fit.log plot.png
