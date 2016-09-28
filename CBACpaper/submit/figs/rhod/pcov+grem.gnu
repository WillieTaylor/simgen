set style line  1 lt rgb "black"
set style line  7 lt rgb "red" 
set style line 31 lt rgb "green"
plot [0:298][0:298] \
'distCB.300' u 3:6 ls 31  pt 7 ps 1.5 not, \
'distCB.300' u 6:3 ls 31  pt 7 ps 1.5 not, \
'psicov.300' u 2:1 ls  7  pt 7 ps 1 not, \
'gremlin.300' u 1:2 ls 7  pt 7 ps 1 not, \
x ls 1 not
set terminal postscript color size 8,7
set output 'plot.ps'
replot
