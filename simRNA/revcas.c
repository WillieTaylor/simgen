/*
cc revcas.c -o revcas util/wt/util.o util/aa/pdbprot.o -lm
*/
#include <stdlib.h>
#include <alloca.h>
#include "util/wt/incl/util.h"
#include "util/wt/incl/geom.h"
#include "util/aa/incl/pdbprot.h"
#include "util/aa/incl/matplot.h"
#include "util/aa/incl/matrix.h"
#define NALLOC 1000
#define NACID 30

main(argc,argv)
int argc; char *argv[];
{
int	i, j, l, n, nout, breaks, anyseq = -1;
Pdbentry_ *prot;
	prot = get_pdb(argv[1],1,1);
	nout = 0;
	for (j=0; j<prot->Chainno; j++)
	{ int	nres = prot->Chains[j].Aano;
		printf("COMPND   %s\n",prot->Compound);
		printf("SOURCE   %s\n",prot->Source);
		for (i=nres-1; i>=0; i--)
		{ Atom_ *Atoms = prot->Chains[j].Atoms;
			nout++;
			printf("ATOM%7d  CA ", nout);
            		printf("%c%3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f",
                		Atoms[i].Alt, aa_code13(Atoms[i].Aa), nout,
                		Atoms[i].X,Atoms[i].Y,Atoms[i].Z,
                		Atoms[i].Occu,Atoms[i].Bfact);
			printf("\n");
		}
	}
	printf("TER\n");
}
