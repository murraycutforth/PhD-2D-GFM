set xrange [0:7]
set yrange [0:3]
set contour base
set cntrparam levels 15
unset surface
set table "pcontTSTM1.dat"
splot "./data/MGFM-CLSVOF-TSTM-MUSCL-211-90-417-prims.dat" u 1:2:7

set xrange [0:7]
set yrange [0:3]
set contour base
set cntrparam levels 15
unset surface
set table "pcontTSTM2.dat"
splot "./data/MGFM-CLSVOF-TSTM-MUSCL-211-90-845-prims.dat" u 1:2:7

set xrange [0:7]
set yrange [0:3]
set contour base
set cntrparam levels 15
unset surface
set table "pcontTSTM3.dat"
splot "./data/VOFGFM-EMOF2-TSTM-MUSCL-211-90-417-prims.dat" u 1:2:7

set xrange [0:7]
set yrange [0:3]
set contour base
set cntrparam levels 15
unset surface
set table "pcontTSTM4.dat"
splot "./data/VOFGFM-EMOF2-TSTM-MUSCL-211-90-845-prims.dat" u 1:2:7
