set term pngcairo
set output 'xprof.png'
set xl 'x'
set yl 'u'
set title 'Vertical profile of x-component of the velocity'
plot[-0.5:0.5] 'xprof.ghia' w p t 'Ghia et al.', \
     'xprof' w l t 'Code'

set output 'yprof.png'
set xl 'y'
set yl 'v'
set title 'Horizontal profile of y-component of the velocity'
plot[-0.5:0.5] 'yprof.ghia' w p t 'Ghia et al.', \
     'yprof' w l t 'Code'

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
