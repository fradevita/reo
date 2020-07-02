set xr [0:pi]
set yr [0:pi]
unset tics
set size ratio 1
set contour
unset surface
set view map
set cntrparam level discr 0.5
set key bottom right
set cntrlabel onecolor
set terminal pngcairo
set output 'plot.png'
splot 'fort.10' i 0 w l lc rgb 'red' t 't = 0', \
      'fort.11' i 0 w l lc rgb 'red' t '', \
      'fort.12' i 0 w l lc rgb 'red' t '', \
      'fort.13' i 0 w l lc rgb 'red' t '', \
      'fort.10' i 1 w l lc rgb 'blue' t 't = T/2', \
      'fort.11' i 1 w l lc rgb 'blue' t '', \
      'fort.12' i 1 w l lc rgb 'blue' t '', \
      'fort.13' i 1 w l lc rgb 'blue' t '', \
      'fort.10' i 2 w l lc rgb 'green' t 't = T', \
      'fort.11' i 2 w l lc rgb 'green' t '', \
      'fort.12' i 2 w l lc rgb 'green' t '', \
      'fort.13' i 2 w l lc rgb 'green' t ''
