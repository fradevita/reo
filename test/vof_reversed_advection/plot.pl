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
splot 'vof' i 0 w l t 't = 0', \
         '' i 1 w l t 't = T/2', \
         '' i 2 w l t 't = T'
