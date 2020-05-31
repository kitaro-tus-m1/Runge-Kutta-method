set terminal png
set output "rk-step2.png"
plot "rk-step2.dat" using 1:2 w l
set output "rk-step2.png"
replot "rk-step2.dat" using 1:3 w l
