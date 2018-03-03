# Create comparison plot of He bubble mass between different GFMs

set terminal epslatex color solid size 16cm,8cm
set output "He_test_3_masserror.tex"

set xlabel 'Time'
set ylabel 'Relative error in bubble mass'

set border lw 3
set key spacing 1.15
set xtics nomirror
set key top left
set ytics 0.05

set xrange [0:280]

plot "./data/OGFM-CLSVOF-shocked_helium_bubble-MUSCL-325-89-masschange.dat"  w l lw 2 title "OGFM", "./data/riemannGFM-CLSVOF-shocked_helium_bubble-MUSCL-325-89-masschange.dat"  w l lw 2 title "RGFM", "./data/MGFM-CLSVOF-shocked_helium_bubble-MUSCL-325-89-masschange.dat"  w l lw 2 title " MGFM"
