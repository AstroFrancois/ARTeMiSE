set xrange [0:2.5] 
set size ratio -1
plot 'q' u 1:2                  
replot 'intersections.dat' u 1:2  
replot 'intersections.dat' u 3:4
replot 'vertex.dat' u 1:2       
