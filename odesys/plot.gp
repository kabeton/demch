set grid x y
show grid
plot "trail.dat" using 1:2 title "y1 computed", "trail.dat" using 1:3 title "y2 computed", "exact.dat" using 1:2 with lines title "y1 from analitical solution", "exact.dat" using 1:3 with lines title "y2 from analitical solution"
pause mouse close
