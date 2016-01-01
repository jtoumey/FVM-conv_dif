reset

set xrange [0:1]
set yrange [0:1]

splot 'phi_distr.dat'
set pm3d map
set samples 100
set isosamples 100
replot
