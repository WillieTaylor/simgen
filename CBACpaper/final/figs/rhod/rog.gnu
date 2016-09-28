set style line 1 lt rgb "black"
set style line 2 lt rgb "red"
set style line 3 lt rgb "green" 
set style line 4 lt rgb "blue"
plot [2:20] \
'films-100-05-8.sort'   u 5:10  ls 2  pt 7 ps 0.5 not, \
'rosetta-100-05-8.sort' u 5:10  ls 3  pt 7 ps 0.5 not, \
'model-100-05-8.sort'   u 5:10  ls 4  pt 7 ps 0.5 not, \
13.9 ls 1 not
set terminal postscript color size 8,7
set output 'plot.ps'
replot
