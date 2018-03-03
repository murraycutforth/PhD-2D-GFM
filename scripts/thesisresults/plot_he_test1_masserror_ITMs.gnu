# Create comparison plot of He bubble mass between different ITMs

set terminal epslatex color solid size 16cm,8cm
set output "He_test_1_masserror.tex"

set xlabel 'Time'
set ylabel 'Relative error in bubble mass'

set border lw 3
set key spacing 1.15
set xtics nomirror
set key bottom left

set xrange [0:280]

plot "./data/OGFM-LS1-shocked_helium_bubble-MUSCL-325-89-masschange.dat"  w l lw 2 title "LS1", "./data/OGFM-LS3-shocked_helium_bubble-MUSCL-325-89-masschange.dat"  w l lw 2 title "LS3", "./data/OGFM-LS5-shocked_helium_bubble-MUSCL-325-89-masschange.dat"  w l lw 2 title "LS5", "./data/OGFM-LS9-shocked_helium_bubble-MUSCL-325-89-masschange.dat"  w l lw 2 title "LS9", "./data/OGFM-CLSVOF-shocked_helium_bubble-MUSCL-325-89-masschange.dat"  w l lw 2 title "CLSVOF" 
