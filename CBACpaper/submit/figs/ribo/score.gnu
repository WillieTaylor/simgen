set style line 1 lt rgb "black"
set style line 2 lt rgb "red"
set style line 3 lt rgb "blue"
set style line 10 lt rgb "green" 
plot [6:14][57:] \
'fullBrms50.dat'   u 2:3  ls 2  pt 7 ps 0.5 not, \
'fullArms50.dat'   u 2:3  ls 3  pt 7 ps 0.5 not
set terminal postscript
set output 'plot.ps'
replot

