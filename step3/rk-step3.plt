set logscale x
set xrange [0.1:1000000]
set yrange [-0.05:1.05]
set xlabel "1+z (= Exp(-N))"
set ylabel "Omega"
plot "rk-step3.dat" using 1:2 w l
replot "rk-step3.dat" using 1:3 w l
replot "rk-step3.dat" using 1:4 w l
set parametric
set trange [-0.05:1.05]
c1 = 1
c2 = 3001
replot c1, t
replot c2, t
