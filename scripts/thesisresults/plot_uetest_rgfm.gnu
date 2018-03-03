# Plot underwater explosion test case for RGFM

set terminal epslatex color solid size 16cm,8cm
set output "UE_RGFM.tex"
set multiplot layout 1,2

set size square
unset xtics
unset ytics
unset key
set border lw 3
unset xlabel
unset ylabel

set xrange [-5:5]
set yrange [-5:5]

set palette rgbformulae 33,13,10

set xlabel '$t = 0.001$'
plot "./data/riemannGFM-LS9-underwater_explosion-MUSCL-400-400-260-prims.dat" w image, "pcontR1.dat" w l lw 2 lc 0

set xlabel '$t = 0.003$'
plot "./data/riemannGFM-LS9-underwater_explosion-MUSCL-400-400-747-prims.dat" w image, "pcontR2.dat" w l lw 2 lc 0
