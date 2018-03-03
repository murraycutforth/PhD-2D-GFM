# Create comparison plot between GFMs on RMI

set terminal epslatex color solid size 16cm,14cm
set output "RMI_comp2.tex"
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

set title 'MGFM with CLSVOF'
plot "./data/MGFM-CLSVOF-RMI-MUSCL-1600-400-6918-prims.dat" u 1:2:3 w image

set title 'VOF-GFM with PY'
plot "./data/VOFGFM-PYVOF-RMI-MUSCL-1600-400-6891-prims.dat" u 1:2:3 w image

set title 'VOF-GFM with EMOF'
plot "./data/VOFGFM-EMOF2-RMI-MUSCL-1600-400-6925-prims.dat" u 1:2:3 w image


