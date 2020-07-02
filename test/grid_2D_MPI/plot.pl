set terminal png
set output 'plot.png'
set title 'boxes'
set key top left
plot 'box1' w p pt 4, 'box2' w p pt 6
