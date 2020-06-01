set terminal postscript enhanced color
set output "rk-step3.eps
set logscale x
set format x "10^{%L}"
set ytics 0.1
set mytics 4
set xrange [0.1:1000000]
set yrange [-0.05:1.05]
set xlabel "1+z (= Exp(-N))"
set ylabel "Omega"
set key invert
set key outside
set key box
set size square
set parametric
set trange [-0.05:1.05]
c1 = 1.0
c2 = 3001.0
c3 = 0.68
plot c1, t lc "dark-gray" notitle
set output "rk-step3.eps"
replot c2, t lc "dark-gray" notitle
set output "rk-step3.eps"
replot t, c3 lc "dark-gray" notitle
set output "rk-step3.eps"
replot "rk-step3.dat" using 1:2 w l lt rgb "blue" lw 2 title "Omega_r"
set output "rk-step3.eps"
replot "rk-step3.dat" using 1:4 w l lt rgb "dark-green" lw 2 title "Omega_m"
set output "rk-step3.eps"
replot "rk-step3.dat" using 1:3 w l lt rgb "red" lw 2 title "Omega_{\lambda}"
