# Create plot of error convergence on CE test with VOF GFM

set terminal epslatex color solid size 14cm,9cm
set output "CEtest_VOFGFM.tex"

set size square

logb(x, base) = log(x)/log(base)
set xlabel '$\log_{2}(N_{\mathrm{cells}})$'


set border lw 3
set key spacing 1.15
set tics nomirror

# set ytics ("\$10^{-6}\$"0.000001, "\$10^{-5}\$"0.00001,"\$10^{-4}\$"0.0001,"\$10^{-3}\$"0.001, "\$10^{-2}\$" 0.01,"\$10^{-1}\$"0.1, "\$10^{0}\$" 1, "\$10^{1}\$" 10, "\$10^{2}\$" 100, "\$10^{3}\$" 1000, "\$10^{4}\$" 10000, "\$10^{5}\$" 100000)
set xrange [5:9]
set yrange [-2.5:-1.0]
set xtics 1
set ytics 0.5
set ylabel '$\log_{10}(L_1)$'
plot -1.0*logb(2,10)*x+0.5 w l lw 6 lc rgb '#00008B' title '\footnotesize $\mathcal{O}(N_{\mathrm{cells}}^{-1})$', "./data/CEtesterror-MGFM.dat" u (logb($1, 2)):(log10($2)) w p pt 7 ps 2 title '\footnotesize MGFM', "./data/CEtesterror-VOFGFM.dat" u (logb($1, 2)):(log10($2)) w p pt 7 ps 2 title '\footnotesize VOF-GFM'

