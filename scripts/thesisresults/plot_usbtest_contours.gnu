set xrange [0:12]
set yrange [0:12]

set contour base
set cntrparam levels 15
unset surface
set table "pcontM3.dat"
splot "./data/MGFM-CLSVOF-underwater_shocked_bubble-MUSCL-400-400-925-prims.dat" u 1:2:7

set contour base
set cntrparam levels 15
unset surface
set table "pcontM4.dat"
splot "./data/MGFM-CLSVOF-underwater_shocked_bubble-MUSCL-400-400-1241-prims.dat" u 1:2:7

set contour base
set cntrparam levels 15
unset surface
set table "pcontR3.dat"
splot "./data/riemannGFM-CLSVOF-underwater_shocked_bubble-MUSCL-400-400-924-prims.dat" u 1:2:7

set contour base
set cntrparam levels 15
unset surface
set table "pcontR4.dat"
splot "./data/riemannGFM-CLSVOF-underwater_shocked_bubble-MUSCL-400-400-1241-prims.dat" u 1:2:7
