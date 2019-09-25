set term png
set output 'xprof.png'
set xl 'x'
set yl 'u'
set title 'Vertical profile of x-component of the velocity'
plot[-0.5:0.5] 'xprof.ghia' w p pt 9 t 'Ghia et al.', \
     'xprof' w l lw 2 t 'reo'

set output 'yprof.png'
set xl 'y'
set yl 'v'
set title 'Horizontal profile of y-component of the velocity'
plot[-0.5:0.5] 'yprof.ghia' w p pt 9 t 'Ghia et al.', \
     'yprof' w l lw 2 t 'reo'

set output 'norm.png'
set xl 'x'
set yl 'y'
set title 'Norm of the velocity for stationary regime'
set size ratio 1
set pm3d
set pm3d map interpolate 10,1
unset key
# jet colormap
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1,	\
0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0, 0.75 1 0.4392 0,	\
0.875 0.9333 0 0, 1 0.498 0 0 )

splot 'out' u 1:2:3

reset
set term pngcairo
set output 'vorticity.png'
set title 'Vorticity contour'
set xl 'x'
set yl 'y'
unset xtics
unset ytics    
unset key
set contour
unset surface
set view map
unset clabel
set cntrparam levels discr -3.0,-2.0,-1.0,-0.5,0.0,0.5,1.0,2.0,3.0,4.0,5.0,6.0
set size ratio 1
splot[1:64][1:64] 'out' u 1:2:5 w l lc rgb 'black'
