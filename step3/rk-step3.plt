set logscale x
set xlabel "-ln(1+z) (= N)"
set ylabel "Omega"
plot "rk-step3.dat" using 1:2
replot "rk-step3.dat" using 1:3
