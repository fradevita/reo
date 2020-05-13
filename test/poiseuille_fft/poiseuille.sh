for cycle in 3 4 5 6; do
    ./poiseuille $cycle > $cycle
done

cat 3 4 5 6 > error
rm 3 4 5 6
gnuplot plot.pl
