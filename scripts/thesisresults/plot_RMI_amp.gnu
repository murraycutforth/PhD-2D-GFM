# Create plot of RMI amplitude

set terminal epslatex color solid size 19cm,10cm
set output "RMI_amp.tex"
set multiplot layout 1,2
set size square

set border lw 3
set key spacing 1.15
set key top right
set tics nomirror
set xrange [0:4]
set xtics 1
set xlabel 'Time'

set ylabel 'Amplitude'
set yrange [0.01:0.11]
set ytics 0.02
plot "./data/VOFGFM-EMOF2-RMI_pert-MUSCL-1600-400-pertamp.dat" w l lw 5 title "VOF-GFM method", "./data/impulsivemodel.dat" u 1:($2 + 0.02) w l lw 5 title "Impulsive model"

set ylabel 'Growth rate'
set yrange [-0.01:0.033]
set ytics 0.01
plot "./data/VOFGFM-EMOF2-RMI_pert-MUSCL-1600-400-pertampderiv.dat" w l lw 5 title "VOF-GFM method", 0.02 w l lw 5 title "Impulsive model"



