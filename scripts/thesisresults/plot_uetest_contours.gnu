set xrange [-5:5]
set yrange [-5:5]
set contour base
set cntrparam levels 15
unset surface
set table "pcontM1.dat"
splot "./data/MGFM-LS9-underwater_explosion-MUSCL-400-400-258-prims.dat" u 1:2:7

set xrange [-5:5]
set yrange [-5:5]
set contour base
set cntrparam levels 15
unset surface
set table "pcontM2.dat"
splot "./data/MGFM-LS9-underwater_explosion-MUSCL-400-400-742-prims.dat" u 1:2:7

set xrange [-5:5]
set yrange [-5:5]
set contour base
set cntrparam levels 15
unset surface
set table "pcontR1.dat"
splot "./data/riemannGFM-LS9-underwater_explosion-MUSCL-400-400-260-prims.dat" u 1:2:7

set xrange [-5:5]
set yrange [-5:5]
set contour base
set cntrparam levels 15
unset surface
set table "pcontR2.dat"
splot "./data/riemannGFM-LS9-underwater_explosion-MUSCL-400-400-747-prims.dat" u 1:2:7
