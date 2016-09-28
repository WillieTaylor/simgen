#include "util/wt/incl/util.h"
#define L 256
#define NREF 500
main () {
FILE	*ful, *srt, *bib, *tex;
char	line[L], mark='%', q='"', cite[NREF][20];
int	tag[NREF], place[NREF];
int	id, nref;
int	i, j, k;
	bib = fopen("select.bib","r");
	tex = fopen("select.tex","w");
	nref = 0;
	while (1) {
		if(read_line(bib,line)<0) break;
		if (line[0]=='@' && line[1]=='s') continue;
		if (strchr(line,'@')) { char tags[L]; int pre,year,suf;
			tag[nref] = (char*)malloc(sizeof(char)*5);
			i = strlen(line)-1;
			line[i] = 0; i--;
			if (isdigit(line[i])) strcpy(tags,line+i-1);
				else strcpy(tags,line+i-2);
			if (tags[2]<'a' || tags[2]>'z') suf = 0;
				else suf = (int)(tags[2]-'a')+1;
			if (tags[0]=='0') pre = 1000; else pre = 0;
			tags[2] = (char)0;
			sscanf(tags,"%d", &year);
			year *= 10;
			tag[nref] = -(pre+year+suf);
			strcpy(cite[nref],strchr(line,'{')+1);
			nref++;
		}
	}
	printf("Sorting %d references into select.tex\n", nref);
	sort(0,0,tag,place,nref,1);
	fprintf(tex,"\\documentstyle[a4,12pt]{report}\n");
	fprintf(tex,"\\begin{document}\n");
	for (i=0; i<nref; i++)
	{ int	pi = place[i];
		fprintf(tex,"REF%d is \\cite{%s}.\n", i,cite[pi]);
	}
	fprintf(tex,"\\clearpage\n");
	fprintf(tex,"\\bibliographystyle{unsrt}\n");
	fprintf(tex,"\\bibliography{select}\n");
	fprintf(tex,"\\end{document}");
}
