#gnuplot
set logscale xy
set xlabel 'N'
set ylabel 'max|e|'
set grid
set title 'Poiseuille convergence'
fit a*x+b 'error' u (log($1)):(log($2)) via a,b
set term pngcairo
set output 'plot.png'
set xtics 8,2,64
plot [6:130]'error' u 1:($2) pt 7 t '', \
     exp(b)*x**a t sprintf("%.0f/N^{%4.2f}", exp(b), -a)
