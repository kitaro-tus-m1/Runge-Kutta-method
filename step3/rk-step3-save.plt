set terminal png
set output "rk-step3.png"
set logscale x
set xlabel "1+z (= Exp(-N))"
set ylabel "Omega"
plot "rk-step3.dat" using 1:2 w l
set output "rk-step3.png"
replot "rk-step3.dat" using 1:3 w l
