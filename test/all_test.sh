# Before executing all tests, clean everything
cd poisson_dirichlet
make clean
cd ..

cd poisson_periodic
make clean
cd ..

cd euler_convergence
make clean
cd ..

cd poiseuille
make clean
cd ..

cd lid
make clean
cd ..

# Run the test for the poisson equation with dirichlet boundary conditions
cd poisson_dirichlet
make poisson_dirichlet
echo "Running poisson dirichlet test"
./poisson_dirichlet
make plot
open plot.png
cd ..

# Run the test for the poisson equation with periodic boundary conditinos
cd poisson_periodic
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

# Run the test for poiseuille
cd poiseuille
make poiseuille
echo "Running poiseuille test"
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
