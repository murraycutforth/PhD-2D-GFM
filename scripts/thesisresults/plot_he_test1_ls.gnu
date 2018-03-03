# Create comparison plot between level set contours of different ITMs

set terminal epslatex color solid size 16cm,14cm
set output "He_test_1_ls.tex"
set multiplot layout 3,1

unset xtics
unset ytics
unset key
set border lw 3
unset xlabel
unset ylabel

set xrange [0:325]
set yrange [0:89]

set palette rgbformulae 33,13,10

set title 'First order upwind'
plot "./data/OGFM-LS1-shocked_helium_bubble-MUSCL-325-89-3083-ls.dat" u 1:2:3 w image, "cont1.dat" w l lw 2 lc 0, "zerocont1.dat" w l lw 5 lc 0

set title 'Fifth order HOUC'
plot "./data/OGFM-LS5-shocked_helium_bubble-MUSCL-325-89-3090-ls.dat" u 1:2:3 w image, "cont5.dat" w l lw 2 lc 0, "zerocont5.dat" w l lw 5 lc 0

set title 'Coupled level set -- VOF'
plot "./data/OGFM-CLSVOF-shocked_helium_bubble-MUSCL-325-89-3089-ls.dat" u 1:2:3 w image, "contC.dat" w l lw 2 lc 0, "zerocontC.dat" w l lw 5 lc 0


