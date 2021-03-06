# Create comparison plot between GFMs on shocked SF6

set terminal epslatex color solid size 16cm,20cm
set output "sf6_comp1.tex"
set multiplot layout 3,1

unset xtics
unset ytics
unset key
set border lw 3
unset xlabel
unset ylabel

set xrange [0.06:0.22]
set yrange [0.06:0.13]
set cbrange [1:20]
set logscale cb

set palette rgbformulae 33,13,10

set title 'MGFM with LS'
plot "./data/MGFM-CLSVOF-shocked_SF6-MUSCL-900-400-884-prims.dat" u 1:2:3 w image

set title 'VOF-GFM with PY'
plot "./data/VOFGFM-PYVOF-shocked_SF6-MUSCL-900-400-880-prims.dat" u 1:2:3 w image 

set title 'VOF-GFM with EMOF'
plot "./data/VOFGFM-EMOF2-shocked_SF6-MUSCL-900-400-878-prims.dat" u 1:2:3 w image

