FC = mpif90
FOPTS = -g -Wall -fcheck=bound -fdefault-real-8 -cpp # -O2 -ffpe-trap=invalid,zero,overflow
CODE_DIR = $(REO)

#HYPRE_DIR = 
INCLUDES = -I$(HYPRE_DIR)/include -I$(CODE_DIR) -I$(FFTW3_DIR)/include
LIBS = -L$(HYPRE_DIR)/lib -lHYPRE -lm -L$(FFTW3_DIR)/lib -lfftw3
FFLAGS = $(FOPTS) $(INCLUDES)
FLAGS = $(OPTS) $(LIBS) -lstdc++

#===============================================================================
# TEST FOR POISSON EQUATION
#===============================================================================
poiseuille: poiseuille.o
	$(FC) -o $@ $< grid_2d.o boundary_conditions.o navier_stokes.o hypre_solver.o navier_stokes_pub.o fft_solver.o volume_of_fluid.o multiphase.o $(FLAGS)

poiseuille.o: poiseuille.f08 grid_2d.o navier_stokes.o navier_stokes_pub.o
	$(FC) $(FFLAGS) -c $<
#===============================================================================

#===============================================================================
# Code
#===============================================================================
navier_stokes.o: $(CODE_DIR)/navier_stokes.f08 grid_2d.o boundary_conditions.o hypre_solver.o navier_stokes_pub.o fft_solver.o volume_of_fluid.o multiphase.o
	$(FC) $(FFLAGS) -c $<

volume_of_fluid.o: $(CODE_DIR)/volume_of_fluid.f08 boundary_conditions.o grid_2d.o
	$(FC) $(FFLAGS) -c $<

multiphase.o: $(CODE_DIR)/multiphase.f08 boundary_conditions.o grid_2d.o volume_of_fluid.o hypre_solver.o
	$(FC) $(FFLAGS) -c $<

hypre_solver.o: $(CODE_DIR)/hypre_solver.f08 grid_2d.o boundary_conditions.o
	$(FC) $(FFLAGS) -c $<

fft_solver.o: $(CODE_DIR)/fft_solver.f08 grid_2d.o boundary_conditions.o
	$(FC) $(FFLAGS) -c $<

boundary_conditions.o: $(CODE_DIR)/boundary_conditions.f08 grid_2d.o navier_stokes_pub.o
	$(FC) $(FFLAGS) -c $<

navier_stokes_pub.o: $(CODE_DIR)/navier_stokes_pub.f08
	$(FC) $(FFLAGS) -c $<

grid_2d.o: $(CODE_DIR)/grid_2d.f08
	$(FC) $(FFLAGS) -c $<

plot:
	gnuplot plot.pl

clean:
	rm -r *.o *.mod poiseuille log error fit.log plot.png
