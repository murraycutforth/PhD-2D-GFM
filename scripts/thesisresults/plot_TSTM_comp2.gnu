# Create comparison plot between GFMs on TSTM

set terminal epslatex color solid size 16cm,8cm
set output "TSTM_comp2.tex"
set multiplot layout 2,2

unset xtics
unset ytics
unset key
set border lw 3
unset xlabel
unset ylabel

set xrange [0:7]
set yrange [0:3]

set palette rgbformulae 33,13,10

set cbtics 1


set cblabel '$\rho$'
set title 'MGFM with CLSVOF'
plot "./data/MGFM-CLSVOF-TSTM-MUSCL-211-90-845-prims.dat" u 1:2:3 w image, "pcontTSTM2.dat" w l lw 2 lc 0

set title 'VOFGFM with EMOF'
plot "./data/VOFGFM-EMOF2-TSTM-MUSCL-211-90-845-prims.dat" u 1:2:3 w image, "pcontTSTM4.dat" w l lw 2 lc 0

set cbtics 0.5

set cblabel '$z$'
set title 'MGFM with CLSVOF'
plot "./data/MGFM-CLSVOF-TSTM-MUSCL-211-90-845-z.dat" u 1:2:3 w image 

set title 'VOFGFM with EMOF'
plot "./data/VOFGFM-EMOF2-TSTM-MUSCL-211-90-845-z.dat" u 1:2:3 w image
