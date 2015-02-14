set key left bottom
set yrange [0:1.1]
plot "temp_distr_i.dat" using 1:2 with linespoints title 'Numerical Solution: Coarse', \
     "temp_distr_ii.dat" using 1:3 with linespoints title 'Analytical Solution', \
     "temp_distr_ii.dat" using 1:2 with linespoints title 'Numerical Solution: Fine'