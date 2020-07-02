set xr [0:1]
set yr [0:1]
set size ratio 1
set contour
#unset surface
#set view map
set cntrparam level discr 0.5
set terminal pngcairo
set output 'plot.png'
set grid
set cntrlabel onecolor
splot 'box1' w l t 'box 1', \
      'box2' w l t 'box 2', \
      'box3' w l t 'box 3', \
      'box4' w l t 'box 4', \
      sqrt((x-0.5)**2 + (y-0.5)**2) - 0.3 + 0.5 w l lc rgb 'black' t 'analytical'
