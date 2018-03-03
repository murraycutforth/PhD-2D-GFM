set xrange [0:0.45]
set yrange [0:0.2]

set contour base
set cntrparam levels 20
unset surface

set table "pcont_sf6_1.dat"
splot "./data/MGFM-LS9-shocked_SF6-MUSCL-900-400-887-prims.dat" u 1:2:7

set table "pcont_sf6_2.dat"
splot "./data/VOFGFM-PYVOF-shocked_SF6-MUSCL-900-400-880-prims.dat" u 1:2:7

set table "pcont_sf6_3.dat"
splot "./data/VOFGFM-EMOF2-shocked_SF6-MUSCL-900-400-878-prims.dat" u 1:2:7
