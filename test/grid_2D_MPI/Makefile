FC = mpif90
FOPTS = -g -Wall -O2 -fcheck=bound -fdefault-real-8 -cpp
CODE_DIR = $(REO)

#HYPRE_DIR = 
INCLUDES =
LIBS =
FFLAGS = $(FOPTS) $(INCLUDES)
FLAGS = $(OPTS) $(LIBS) -lstdc++

#===============================================================================
# TEST FOR CHECK GRID
#===============================================================================
check_grid: check_grid.o
	$(FC) -o $@ $< grid_2d.o $(FLAGS)

check_grid.o: check_grid.f08 grid_2d.o
	$(FC) $(FFLAGS) -c $<
#===============================================================================

#===============================================================================
# Code
#===============================================================================
grid_2d.o: $(CODE_DIR)/grid_2d.f08
	$(FC) $(FFLAGS) -c $<

plot:
	gnuplot plot.pl

clean:
	rm -r *.o *.mod check_grid box* plot.png
