set logscale x
set xlabel "1+z (= Exp(N))"
set ylabel "Omega"
plot "rk-step3.dat" using 1:2
replot "rk-step3.dat" using 1:3
