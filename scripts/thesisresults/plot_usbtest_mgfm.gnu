# Plot underwater shocked bubble test case for MGFM

set terminal epslatex color solid size 16cm,8cm
set output "USB_MGFM.tex"
set multiplot layout 1,2

set size square
unset xtics
unset ytics
unset key
set border lw 3
unset xlabel
unset ylabel

set xrange [0:12]
set yrange [0:12]

set palette rgbformulae 33,13,10

set xlabel '$t = 0.0125$'
plot "./data/MGFM-CLSVOF-underwater_shocked_bubble-MUSCL-400-400-925-prims.dat" w image, "pcontR3.dat" w l lw 2 lc 0

set xlabel '$t = 0.0375$'
plot "./data/MGFM-CLSVOF-underwater_shocked_bubble-MUSCL-400-400-1241-prims.dat" w image, "pcontR4.dat" w l lw 2 lc 0
