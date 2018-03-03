# Plot underwater explosion test case for MGFM

set terminal epslatex color solid size 16cm,8cm
set output "UE_MGFM.tex"
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
plot "./data/MGFM-LS9-underwater_explosion-MUSCL-400-400-258-prims.dat" w image, "pcontM1.dat" w l lw 2 lc 0

set xlabel '$t = 0.003$'
plot "./data/MGFM-LS9-underwater_explosion-MUSCL-400-400-742-prims.dat" w image, "pcontM2.dat" w l lw 2 lc 0
