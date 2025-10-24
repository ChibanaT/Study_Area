# plot_mesh.plt - use: gnuplot -persist plots/plot_mesh.plt
set title "Annulus mesh (mm) - ThermoInduce"
set size ratio -1
set xlabel "x (mm)"
set ylabel "y (mm)"
set grid

# draw the ring lines (blocks in mesh.dat)
plot "mesh.dat" index 0 with lines lt 1 lw 1 notitle, \
     "" index 1 with lines lt 1 lw 1 notitle, \
     "" index 2 with lines lt 1 lw 1 notitle, \
     "" index 3 with lines lt 1 lw 1 notitle, \
     "" index 4 with lines lt 1 lw 1 notitle, \
     "" index 5 with lines lt 1 lw 1 notitle, \
     "" index 6 with lines lt 1 lw 1 notitle, \
     "" index 7 with lines lt 1 lw 2 lc rgb 'red' title 'radial spokes'
# If Nr is larger, you can extend indices or use: plot for [i=0:15] 'mesh.dat' index i with lines
