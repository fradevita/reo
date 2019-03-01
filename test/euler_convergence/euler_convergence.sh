for cycle in 3 4 5 6 7 8; do
    ./euler_convergence $cycle > $cycle
done

cat 3 4 5 6 7 8 > error
rm 3 4 5 6 7 8
gnuplot plot.pl
