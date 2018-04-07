set xrange [1:3]
set yrange [0.5:1]
set contour base
set cntrparam levels 20
unset surface
set table "pcontRMI1.dat"
splot "./data/MGFM-CLSVOF-RMI-MUSCL-1600-400-3402-prims.dat" u 1:2:7

set xrange [1:3]
set yrange [0.5:1]
set contour base
set cntrparam levels 20
unset surface
set table "pcontRMI2.dat"
splot "./data/VOFGFM-PYVOF-RMI-MUSCL-1600-400-3401-prims.dat" u 1:2:7

set xrange [1:3]
set yrange [0.5:1]
set contour base
set cntrparam levels 20
unset surface
set table "pcontRMI3.dat"
splot "./data/VOFGFM-EMOF2-RMI-MUSCL-1600-400-3396-prims.dat" u 1:2:7
