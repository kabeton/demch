set grid x y
show grid
plot "trail.dat" using 1:2 title "y1 computed", "trail.dat" using 1:3 title "y2 computed"
pause mouse close
