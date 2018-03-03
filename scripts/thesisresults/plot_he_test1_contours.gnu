set xrange [0:325]
set yrange [0:89]
set contour base
set cntrparam levels discrete 0.0
unset surface
set table "zerocont1.dat"
splot "./data/OGFM-LS1-shocked_helium_bubble-MUSCL-325-89-3083-ls.dat"

set xrange [0:325]
set yrange [0:89]
set contour base
set cntrparam levels auto 25
unset surface
set table "cont1.dat"
splot "./data/OGFM-LS1-shocked_helium_bubble-MUSCL-325-89-3083-ls.dat"

set xrange [0:325]
set yrange [0:89]
set contour base
set cntrparam levels discrete 0.0
unset surface
set table "zerocont5.dat"
splot "./data/OGFM-LS5-shocked_helium_bubble-MUSCL-325-89-3090-ls.dat"

set xrange [0:325]
set yrange [0:89]
set contour base
set cntrparam levels auto 25
unset surface
set table "cont5.dat"
splot "./data/OGFM-LS5-shocked_helium_bubble-MUSCL-325-89-3090-ls.dat"

set xrange [0:325]
set yrange [0:89]
set contour base
set cntrparam levels discrete 0.0
unset surface
set table "zerocontC.dat"
splot "./data/OGFM-CLSVOF-shocked_helium_bubble-MUSCL-325-89-3089-ls.dat"

set xrange [0:325]
set yrange [0:89]
set contour base
set cntrparam levels auto 25
unset surface
set table "contC.dat"
splot "./data/OGFM-CLSVOF-shocked_helium_bubble-MUSCL-325-89-3089-ls.dat"
