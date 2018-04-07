# Create comparison plot between GFMs on RMI

set terminal epslatex color solid size 16cm,14cm
set output "RMI_comp1.tex"
set multiplot layout 3,1

unset xtics
unset ytics
unset key
set border lw 3
unset xlabel
unset ylabel

set xrange [1:3]
set yrange [0.5:1]

set palette rgbformulae 33,13,10

set cblabel '$\rho$'

set title 'MGFM with CLSVOF'
plot "./data/MGFM-CLSVOF-RMI-MUSCL-1600-400-3402-prims.dat" u 1:2:3 w image, "pcontRMI1.dat" w l lw 2 lc 0

set title 'VOF-GFM with PY'
plot "./data/VOFGFM-PYVOF-RMI-MUSCL-1600-400-3401-prims.dat" u 1:2:3 w image, "pcontRMI2.dat" w l lw 2 lc 0

set title 'VOF-GFM with EMOF'
plot "./data/VOFGFM-EMOF2-RMI-MUSCL-1600-400-3396-prims.dat" u 1:2:3 w image, "pcontRMI3.dat" w l lw 2 lc 0

