#include "util.hpp"
#include "geom.hpp"
#include "cell.hpp"
#include "data.hpp"

void paramin ( char*, int );

int params ( FILE *run ) {
// read the global behaviour values and model data
// NB commands must come in the order specified below as there is no loop 
int	n, io, nmodels;
float	s, sget, sput;
char	line[222];
	Data::norun = 0;
	Data::noview = 0;
	Data::putpdb = 0;
	Data::shrink = 0.0;	// used as 1-shrink (-ve = mass factor)
	Data::scalein = 0.1;	// default unless specified by SCALE <in> <out>
	Data::scaleout = 1.0;	// unless specified by SCALE
	LOOP {
		io = read_line(run,line);
		LOOP {
			if (io < 0) { Pt(Unexpected end of file\n) exit(1); }
			if (io==0) continue;
			if (line[0]=='#') {		// skip comment
				io = read_line(run,line);
				continue;
			}
			break;
		}
		if (line[0]=='N' && line[2]=='R') {	// NORUN flag (skip driver)
			Data::norun = 1;
			continue;
		}
		if (line[0]=='N' && line[2]=='M') {	// NOMOVE flag (no linker,shaker,bumper 
			Data::norun = 2;			//	or driver)
			continue;
		}
		if (line[0]=='N' && line[2]=='V') {	// NOVIEW flag (skip viewer)
			Data::noview = 1;
			continue;
		}
		if (line[0]=='P' && line[2]=='T') {	// PUTPDB flag (dump pdb on exit)
			Data::putpdb = 1;
			continue;
		}
		if (line[0]=='H' && line[2]=='D') {	// HIDDEN [value] (hide object contents)
			Data::hidden = 1;
			if (io > 7) {
				sscanf(line+6,"%d", &n);
				Data::hidden = n;
			}
			continue;
		}
		if (line[0]=='T' && line[2]=='N') {	// TINKER [update] (apply local refinement)
			Data::tinker = -10;		// -10 = refine starting geom ever 10th frame
			if (io > 7) {			// n>0 = refine updated geom every nth cycle
				sscanf(line+6,"%d", &n);	// n<0 keep start geom
				Data::tinker = n;
			}
			continue;
		}
		if (line[0]=='S' && line[1]=='C') { // output SCALE factor (saved in world->far)
			Ps(line) NL
			sscanf(line+6,"%f %f", &sget, &sput);
			Data::scalein = sget;
			Data::scaleout = sput;
			continue;
		}
		if (line[0]=='S' && line[2]=='R') { // SHRINK factor (saved later in scene->far)
			Ps(line) NL
			sscanf(line+7,"%f", &s);
			Data::shrink = s;
			continue;
		}
		break;
	} // no pre-commands found so time to read params
	if (line[0] != 'P') {
		printf("Need: PARAM filename\n");
		return 2;
	}
	nmodels = 0;
	while (1) { char atomtype; // read in param file(s)
		printf("Reading parameters from %s\n", line+6);
		paramin(line+6,nmodels);
		nmodels++;
		read_line(run,line);
		if (line[0] == 'P') continue;
		if (line[0] == 'M' || line[0] == 'E') { // MODEL or END to finish
			printf("End of parameters\n\n");
			Pt(Read) Pi(nmodels) NL
			DO(i,nmodels) // set the atomic bond lengths
			{ float len = Data::model[i].sizes[Data::depth]
				    + Data::model[i].bonds[Data::depth];
				printf("Model %d is", i);
				if (Data::model[i].moltype==0) { // protein
					Data::bondCA = len;
					Data::scalein = len/3.8;
					Data::scaleout = 3.8/len;
					Pt(protein with) Pr(Data::bondCA) Pr(Data::scalein) NL
				}
				if (Data::model[i].moltype==1) { // nucleic
					Data::bondPP = len;
					Pt(nucleic with) Pr(Data::bondPP) NL
				}
				if (Data::model[i].moltype==2) { // chemistry
					Pt(chemical) NL
				}
				if (Data::model[i].moltype==3) { // cells
					Pt(cells) NL
				}
			}
			break;
		}
		Pt(Bad) Ps(line) NL exit(1);
	}
	Data::nmodels = nmodels;
	if (line[0] == 'M') { int model;
		sscanf(line+6,"%d", &model);
		printf("Using model %d\n", model);
	}
}

void paramin ( char *param, int n ) {
int	i, j, io;
FILE	*dat;
char	*at;
char	line[222];
char	name;
Data	*model = Data::model+n;
	dat = fopen(param,"r");
	if (!dat) { printf("%s parameter file not found\n", param); exit(1); }
	io = read_line(dat,line);
	model->moltype = model->subtype = 0; // default (any protein-like chain)
        if (line[0]=='P') model->subtype = 1; // PROTein (CA-CA = 3.8)
        if (line[0]=='R') { model->moltype = 1; model->subtype = 0; } // RNA
        if (line[0]=='D') { model->moltype = 1; model->subtype = 1; } // DNA
        if (line[1]=='H') model->moltype = 2; // CHEM (needs fixing)
        if (line[1]=='E') model->moltype = 3; // CELLs
	model->colours = 0;		// default = colour by level
	if (strlen(line) > 4) {
        	if (toupper(line[5])=='S') model->colours = 1; // PROT with SSE coloured red/green
        	if (toupper(line[5])=='B') model->colours = 2; // flash Bumps in a different colour
	}
	name = line[io-1];
	i = 0;
	while (1)
	{ int	io = read_line(dat,line);
		if (io < 0) break;
		if (io == 0) continue;
		if (line[0] == '/') continue;
		*(strstr(line,"/")) = '\0';
		at = line;
		for (j=0; j<N; j++) { int in;
			sscanf(at,"%d", &in);
			if (in != 0) Data::depth = max(Data::depth,j);
			printf("%5d ", in);
			at = strstr(at,",");
			if (!at) break;
			at++;
			switch (i) {
				case 0 :
					model->shape[j] = in;
					//	-ve = render as wireframe
					break;
				case 1 :
					model->sizes[j] = 0.1*(float)abs(in);
					//	-ve = use values in prox and dist 
					if (in<0) model->local[j] = 1; else  model->local[j] = 0;
					break;
				case 2 :
					model->bumps[j] = 0.1*(float)in;
					//	-ve
					break;
				case 3 :
					model->links[j] = in;
					//	-ve
					break;
				case 4 :
					model->chain[j] = in;
					//	-ve = circular
					break;
				case 5 :
					model->kicks[j] = 0.01*(float)in;
					//	-ve
					break;
				case 6 :
					model->keeps[j] = 0.01*(float)in;
					//	-ve = shell for types 1,3 or tube = open
					break;
				case 7 :
					model->bonds[j] = 0.1*(float)in;
					//	-ve = no bond between cousins
					if (in<0) model->split[j] = 1; else  model->split[j] = 0;
					break;
				case 8 :
					model->repel[j] = 0.01*(float)in; // hard bump
					//	-ve don't repel children within parent
					break;
				case 9 :
					model->rejel[j] = 0.001*(float)in; // soft bump
					//	-ve don't repel children within parent
					break;
			}
		} NL
		i++;
		if (i==N) break;
	}
	fclose(dat);
}
/*
main () {
	paramin("xxx");
	Pi(Data::model[0].moltype) NL
	Pi(Data::model[0].shape[3]) NL
}
*/
