set xr [0:1]
set yr [0:1]
set size ratio 1
set contour
unset surface
set view map
set cntrparam level discr 0.5
set terminal pngcairo
set output 'plot.png'
set grid
set cntrlabel onecolor
splot 'h' w l t 'reconstructed', sqrt((x-0.5)**2 + (y-0.5)**2) - 0.433 + 0.5 t 'analytical'
