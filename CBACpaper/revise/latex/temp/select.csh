#! /bin/csh
#
awk -f select.awk $argv[1].bib >! select.bib
./select
cat head.bib select.bib >! temp.bib
mv temp.bib select.bib
latex select
bibtex select
latex select
mv select.bbl select.bbl.tex
latex selected
dvips selected
gv selected
