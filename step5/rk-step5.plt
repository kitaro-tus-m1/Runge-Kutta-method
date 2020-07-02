set logscale x
set format x "10^{%L}"
set ytics 0.1
set mytics 5
set xrange [0.1:1000000]
set yrange [-0.05:1.05]
set xlabel "1+z  ( = e^{-N} )"
set ylabel "{/Symbol W}"
set key outside invert reverse box Left width 2
set size square
set parametric
set trange [-0.05:1.05]
c1 = 1.0
c2 = 3001.0
c3 = 0.68
plot c1, t lc "dark-gray" notitle
replot c2, t lc "dark-gray" notitle
replot t, c3 lc "dark-gray" notitle
replot "rk-step5.dat" using 1:2 w l lt rgb "blue" title "{/Symbol W}_r"
replot "rk-step5.dat" using 1:4 w l lt rgb "dark-green" title "{/Symbol W}_m"
replot "rk-step5.dat" using 1:3 w l lt rgb "red" title "{/Symbol W}_{/Symbol L}"
