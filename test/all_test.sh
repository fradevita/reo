# Run the test for the poisson equation with dirichlet boundary conditions
cd poisson_dirichlet
make poisson_dirichlet
echo "Running poisson dirichlet test"
./poisson_dirichlet
make plot
open plot.png
cd ..

# Run the test for the poisson equation with periodic boundary conditinos
cd poisson_periodic_mgr
make poisson_periodic
echo "Running poisson periodic test"
./poisson_periodic
make plot
open plot.png
cd ..

# Run the test for the convergence of euler equation
cd euler_convergence
make euler_convergence
echo "Running euler convergence test"
sh euler_convergence.sh
open plot.png
cd ..

# Run the test for poiseuille periodic with FFT
cd poiseuille_fft
make poiseuille
echo "Running poiseuille fft test"
sh poiseuille.sh
open plot.png
cd ..

# Run the test for poiseuille periodic with HYPRE
cd poiseuille_mgr
make poiseuille
echo "Running poiseuille multigrid test"
sh poiseuille.sh
open plot.png
cd ..

# Run the test for poiseuille inflow / outflow
cd poiseuille_io
make poiseuille
echo "Running poiseuille inflow multigrid test"
sh poiseuille.sh
open plot.png
cd ..

# Run the driven cavity test
cd lid
make lid
echo "Running lid driven cavity test"
./lid
make plot
open *.png
cd ..

echo "DONE"
