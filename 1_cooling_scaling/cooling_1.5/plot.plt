set ytics
set y2tics
set grid

set xrange [0:2e-10]
set yrange [-1.2:1.2]
set y2range [0:900]

set xlabel "Time (s)"
set ylabel "Magnetisation"
set y2label "Temperature (K)"

plot "output" u 1:8 w l title "Mag_{net}", "output" u 1:($15*$16) w l title "Mag_2", "output" u 1:($11*$12) w l title "Mag_1", "output" u 1:2 w l title "T" axes x1y2, "output" u 1:3 title "T_L" axes x1y2 w l, "output" u 1:4 w l title "T_S" axes x1y2,

pause -1
