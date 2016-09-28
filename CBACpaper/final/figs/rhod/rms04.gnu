set style line 1 lt rgb "black"
set style line 2 lt rgb "red"
set style line 3 lt rgb "blue"
set style line 10 lt rgb "green" 
plot [2:10][130:] \
'films-100-04-8.sort'   u 5:15  ls 2  pt 7 ps 0.5 not, \
'model-100-04-8.sort'   u 5:15  ls 3  pt 7 ps 0.5 not, \
186.25 ls 10 not
set terminal postscript color size 8,7
set output 'plot.ps'
replot
