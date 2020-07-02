#gnuplot
set pm3d 
set pm3d map interpolate 10,2
set size ratio -1
set terminal pngcairo
set output 'plot.png'
unset key
splot 'box1', 'box2', 'box3'
