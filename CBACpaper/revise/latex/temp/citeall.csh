#! /bin/csh
#
rm current.bbl
latex current
bibtex current
latex current
mv current.bbl current.bbl.tex
latex sorted
dvips sorted
