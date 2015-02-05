set key left bottom
set yrange [0:1.1]
plot "temp_distr.dat" using 1:2 with line title 'Numerical Solution', \
     "temp_distr.dat" using 1:3 with linespoints title 'Analytical Solution'