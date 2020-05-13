#gnuplot
set logscale xy
set xlabel 'N'
set ylabel 'max|e|'
set grid
set title 'Poisson periodic multigrid'
fit a*x+b 'error' u (log($1)):(log($2)) via a,b
set term png
set output 'plot.png'
set xtics 8,2,256
plot [8:128]'error' u 1:($2) pt 7 t '', \
     exp(b)*x**a t sprintf("%.0f/N^{%4.2f}", exp(b), -a)
