reset
set pm3d map
set colorsequence classic
set contour base
unset key
set autoscale fix
set cntrparam levels discr 0.0
unset clabel
set isosamples 200,200
set palette model RGB
set palette defined
set palette rgb 3,13,10
set cbrange [-4:4]
splot 'oooo' i 1480 u 1:2:6  with pm3d nocontour, \
       'fort.101' with lines lc rgb 'black' nosurface



