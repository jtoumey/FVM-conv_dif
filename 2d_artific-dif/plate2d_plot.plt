set pm3d map
set cbrange [0:110]
set palette defined (0 "blue",17 "#00ffff",33 "white",50 "yellow",\
    66 "red",100 "#990000",101 "grey")
splot 'phi_distr.dat' matrix
