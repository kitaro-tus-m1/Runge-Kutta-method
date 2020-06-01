set logscale x
set format x "10^{%L}"
set ytics 0.1
set mytics 4
set xrange [0.1:1000000]
set yrange [-0.05:1.05]
set xlabel "1+z (= Exp(-N))"
set ylabel "Omega"
plot "rk-step3.dat" using 1:2 w l
replot "rk-step3.dat" using 1:3 w l
replot "rk-step3.dat" using 1:4 w l
set arrow from 0.1,0.68 to 1,0.68 nohead
set parametric
set trange [-0.05:1.05]
c1 = 1.0
c2 = 3001.0
replot c1, t
replot c2, t
