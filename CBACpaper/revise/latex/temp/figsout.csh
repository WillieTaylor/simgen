sed '/begin/,$ d' paper.tex > figsout.tex
echo "\\begin{document}" >> figsout.tex
echo "" >> figsout.tex
echo "\\setcounter{figure}{0}" >> figsout.tex
echo "" >> figsout.tex
foreach file (`cat figsout.list`)
	echo $file
	echo "" >> figsout.tex
	echo "%" $file.tex >> figsout.tex
	echo "" >> figsout.tex
	grep -v '^\%' $file.tex > temp.tex
	awk '{if(match($1,"begin{figure")){p=1;}; if(match($1,"caption")){printf("\\caption{}\n\\end{figure}\n\\clearpage\n\\setcounter{subfigure}{0}\n"); p=0;}; if(p==1){print $0}}' temp.tex >> figsout.tex
	echo "" >> figsout.tex
end
echo "\\end{document}" >> figsout.tex
rm temp.tex
