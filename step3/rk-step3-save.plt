set terminal png
set output "rk-step3.png"
set logscale x
set xrange [0.1:1000000]
set yrange [-0.05:1.05]
set xlabel "1+z (= Exp(-N))"
set ylabel "Omega"
plot "rk-step3.dat" using 1:2 w l
set output "rk-step3.png"
replot "rk-step3.dat" using 1:3 w l
set output "rk-step3.png"
replot "rk-step3.dat" using 1:4 w l
set parametric
set trange [-0.05:1.05]
const = 3001
set output "rk-step3.png"
replot const, t
