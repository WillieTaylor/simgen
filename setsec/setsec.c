/*
cc -c lines.c
cc setsec.c -o setsec lines.o util/gradproj.o util/aa/stutest.o util/aa/bestrot.o util/wt/util.o util/wt/geom.o util/wt/sort.o util/aa/pdbprot.o util/aa/matrix.o util/aa/siva.o util/aa/ql.o -lm

setsec 3chy.pdb 1.3 1.0 0.25
rasmol -script test.ras test.out

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

typedef struct {
        int     a, b;
        float   c;
} Pairs;

typedef struct {
	char	pdb, rms, seg;
} Secs;

float cadist();
float dmatch();
int   axisfit();

main(argc,argv) int argc; char *argv[]; {
Vec	cas[999], axes[5];
Pairs	betaij[999];
float	gapen, betaw, betae;
int	i, j, l, n, nout, breaks, anyseq = -1;
char	code[20], line[100], secs[10000], pdbsec[10000], chid = ' ';
/* NB 10000 used for v.big res numbers */
Secs	*secstr;
Pdbentry_ *prot;
FILE	*secout, *pdbout, *simout;
int	fit;
	sscanf(argv[2],"%f", &gapen);
	sscanf(argv[3],"%f", &betaw);
	sscanf(argv[4],"%f", &betae);
	Pr(gapen) Pr(betaw) Pr(betae) NL
	secout = fopen("secs.dat","w");
/*
	if (argv[1][0]=='>') strcpy(code,argv[1]+1);
		else	     strcpy(code,argv[1]);
	l = strlen(code);
        if (l==5) { chid = code[l-1]; code[l-1] = 0; }
	chid = toupper(chid);
*/
//        strcpy(line,"pdb/");
//        strcat(line,argv[1]);
        strcpy(line,argv[1]);
	prot = get_pdb(line,1,1);
	code[4] = 0;
	nout = 0;
	for (j=0; j<prot->Chainno; j++)
	{ int	nres = prot->Chains[j].Aano,
		nsec = prot->Chains[j].Secsno,
		inseg[99], nseg, in;
	  float asum, bsum, csum, apct, bpct, cpct,
		s,t, x,y,z;
		if (nres==0) continue;
/*
		if (prot->Chains[j].Chid != chid) continue;
*/
		if ((anyseq>-1) && (prot->Chains[j].Chid == chid)) /* the next chain continues previous */
		{ int	difs = 0;			/* unless it is identical */
			for (i=0; i<min(30,nres); i++) { /* check first 30 res */
				if (prot->Chains[j].Atoms[i].Aa == prot->Chains[anyseq].Atoms[i].Aa) 
					continue;
				difs = 1;
				break;
			}
			if (difs) {
				for (i=0; i<nres; i++)
				{ Atom_ *Atoms = prot->Chains[j].Atoms;
					nout++;
					printf("ATOM%7d  CA ", nout);
            				printf("%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f",
                				Atoms[i].Alt, aa_code13(Atoms[i].Aa),
                				chid,Atoms[i].Resno,Atoms[i].Rid,
                				Atoms[i].X,Atoms[i].Y,Atoms[i].Z,
                				Atoms[i].Occu,Atoms[i].Bfact);
					printf("\n");
				}
				continue;
			}
		}
		if (anyseq>-1) break;
		n = 0;
		line[n] = ' ';
		for (i=0; i<strlen(prot->Compound); i++) {
			if (n==99) break;
			if (isspace(prot->Compound[i]) && isspace(line[n])) continue;
			n++;
			line[n] = prot->Compound[i];
		}
		line[n] = 0;
		pdbout = fopen("secout.pdb","w");
		fprintf(pdbout,"COMPND   %s\n",line);
		n = 0;
		line[n] = ' ';
		for (i=0; i<strlen(prot->Source); i++) {
			if (n==99) break;
			if (isspace(prot->Source[i]) && isspace(line[n])) continue;
			n++;
			line[n] = prot->Source[i];
		}
		line[n] = 0;
		fprintf(pdbout,"SOURCE   %s\n",line);
		fprintf(pdbout,"REMARK    ");
		chid = prot->Chains[j].Chid;
		fprintf(pdbout,"%s %c",prot->Pdbcode, chid);
		fprintf(pdbout,"%6.2f ",prot->Resol);
		fprintf(pdbout," %d-%d ",  prot->Chains[j].Atoms[0].Resno,
				   	prot->Chains[j].Atoms[nres-1].Resno);
		fprintf(pdbout," %d ",nres);
		fprintf(pdbout,"  PERCENT SEC.STR ");
		secstr = (Secs*)alloca(sizeof(Secs)*(nres+2));
		for (i=0; i<10000; i++) secs[i] = pdbsec[i] = 'x';
		for (i=0; i<=nres; i++) secstr[i].pdb = secstr[i].rms = secstr[i].seg = '-';
		asum = bsum = csum = -1.0;
		if (nsec) {
			for (i=0; i<2000; i++) secs[i] = pdbsec[i] = '-';
			asum = bsum = csum = 0.0;
			for (i=0; i<nsec; i++)
			{ Secstr_ *ss = prot->Chains[j].Secs+i;
			  int	len = ss->End - ss->Beg + 1;
			  char  t='-'; int k;
				if (ss->Sectype == HELIX) {
					if (ss->Type==1) { t = 'A'; asum += (float)len; }
					if (ss->Type==5) { t = '3'; csum += (float)len; }
				}
				if (ss->Sectype == SHEET) { t = 'B'; bsum += (float)len; }
				if (ss->Beg < 0) ss->Beg = 0;
				if (ss->End < 0) ss->Beg = 1;
				for (k=ss->Beg; k<=ss->End; k++) pdbsec[k] = t;
				if (t=='A') pdbsec[ss->Beg] = pdbsec[ss->End] = 'a';
				if (t=='B') pdbsec[ss->Beg] = pdbsec[ss->End] = 'b';
			}
			asum = 100.0*asum/(float)nres;
			bsum = 100.0*bsum/(float)nres;
			csum = 100.0*csum/(float)nres;
		}
		if (asum>99.9) asum = bsum = -1.0;
		if (bsum>99.9) asum = bsum = -1.0;
		printf("%5.1f%5.1f ", asum,bsum);
		breaks = consec_seen(prot->Chains+j,secs,&apct,&bpct,&cpct,secstr);
		fprintf(pdbout,"(%5.1f%5.1f) ", apct,bpct);
		fprintf(pdbout,"\n");
		lineup(argv[1],prot,j,secstr,gapen,betaw,betae,betaij);
		bsum = 0.0;
		for (i=0; i<nres; i++) bsum += prot->Chains[j].Atoms[i].Bfact;
		if (secstr[1].seg=='-' && secstr[2].seg!='-') secstr[1].seg=secstr[2].seg;
		if (secstr[nres].seg=='-' && secstr[nres-1].seg!='-') secstr[nres].seg=secstr[nres-1].seg;
		for (i=1; i<=nres; i++) printf("%c", secstr[i].seg); NL
		fprintf(pdbout,"REMARK    Length =%4d\n", nres);
		t = 0.0;
		nseg = in = 0;
		for (i=0; i<nres; i++)
		{ Atom_ *Atoms = prot->Chains[j].Atoms;
		  char c = toupper(secstr[i+1].seg);
			s = 1.0;
			if (c=='A') s = 3.0;
			if (c=='3') s = 2.5;
			if (c=='B') s = 2.0;
			if (t != s && t > 0.0) {
				inseg[nseg] = in;
				nseg++;
				in = 0;
			}
			t = Atoms[i].Bfact;
			/* t = 2.0*exp(-t*t*0.001)-1.0; */
			nout++;
			fprintf(pdbout,"ATOM%7d  CA ", nout);
            		fprintf(pdbout," %3s A%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f",
                		aa_code13(Atoms[i].Aa),
                		Atoms[i].Resno,Atoms[i].Rid,
                		Atoms[i].X,Atoms[i].Y,Atoms[i].Z,
                		t,s);
			fprintf(pdbout,"\n");
			t = s;
			in++;
			cas[i+1].x=Atoms[i].X; cas[i+1].y=Atoms[i].Y; cas[i+1].z=Atoms[i].Z;
			// cas runs 1..N
		}
		inseg[nseg] = in;
		anyseq = j;
		fprintf(secout,">%s  %d %6.2f\n", argv[1],nres,prot->Resol);
		for (i=0; i<nres; i++)
		{ Atom_ *Atoms = prot->Chains[j].Atoms;
		  int	at = Atoms[i].Resno;
			if (at < 0) at = 0;
			secstr[i+1].pdb = pdbsec[at];
			fprintf(secout,"%c", Atoms[i].Aa);
		} 
		fprintf(secout,"\n");
		for (i=1; i<=nres; i++) fprintf(secout,"%c", secstr[i].pdb); fprintf(secout,"\n");
		for (i=1; i<=nres; i++) fprintf(secout,"%c", secstr[i].rms); fprintf(secout,"\n");
		for (i=1; i<=nres; i++) fprintf(secout,"%c", secstr[i].seg); fprintf(secout,"\n");
		fclose(secout);
/*
		n = 0;
		for (i=1; i<=nres; i++) if (secstr[i].pdb==secstr[i].rms) n++;
		printf("pdb/rms = %5.1f\n", 100.0*(float)n/(float)nres);
		n = 0;
		for (i=1; i<=nres; i++) if (secstr[i].pdb==secstr[i].seg) n++;
		printf("pdb/seg = %5.1f\n", 100.0*(float)n/(float)nres);
		n = 0;
		for (i=1; i<=nres; i++) if (secstr[i].seg==secstr[i].rms) n++;
		printf("seg/rms = %5.1f\n", 100.0*(float)n/(float)nres);
*/
		simout = fopen("secout.dat","w");
		fit = axisfit(cas,1,nres,axes);
		Pv(axes[0]) NL // CoG
		Pv(axes[1]) NL // A-axis
		Pv(axes[2]) NL // B-axis
		Pv(axes[3]) NL // C-axis
		Pv(axes[4]) NL // axis lengths
		s = 10.0;
		if (fit > 10) { // oblate
			x = axes[1].x*s; y = axes[1].y*s; z = axes[1].z*s;
		} else {	// prolate
			x = axes[3].x*s; y = axes[3].y*s; z = axes[3].z*s;
		}
		fprintf(simout,"GROUP %d %d   %4.1f %4.1f %4.1f\n", fit, nseg+1, x,y,z);
		t = -1.0;
		nseg = 0;
		nout = 0;
		for (i=0; i<nres; i++)
		{ Atom_ *Atoms = prot->Chains[j].Atoms;
		  char c = toupper(secstr[i+1].seg);
			s = 1.0;
			if (c=='A') s = 3.0;
			if (c=='3') s = 2.5;
			if (c=='B') s = 2.0;
			if (t != s) { int sort, ins = inseg[nseg];
				if (s > 2.5) sort = 1;
				if (s < 2.5) sort = 2;
				if (s < 1.5) sort = 0;
				if ( sort == 0) {
					fprintf(simout,"\tGROUP %d %d\n", sort, ins);
				} else {
					fit = axisfit(cas,i+1,i+ins,axes);
					t = 10.0;
					x = axes[1].x*t; y = axes[1].y*t; z = axes[1].z*t;
					fprintf(simout,"\tGROUP %d %d   %4.1f %4.1f %4.1f\n", sort, ins, x,y,z);
				}
				nseg++;
			}
			t = 1.0;
			nout++;
			fprintf(simout,"\t\tATOM%7d  CA ", nout);
            		fprintf(simout," %3s A%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f",
                		aa_code13(Atoms[i].Aa),
                		Atoms[i].Resno,Atoms[i].Rid,
                		Atoms[i].X,Atoms[i].Y,Atoms[i].Z,
                		t,s);
			fprintf(simout,"\n");
			t = s;
		}
		fprintf(simout,"SHEET pdbid\n");
		for (i=0; i<999; i++) {
			if (betaij[i].c < 0.0) break;
			fprintf(simout,"BETA %5d %5d %5d\n", betaij[i].a, betaij[i].b, (int)(betaij[i].c+NOISE));
		}
	}
	if (anyseq>-1) fprintf(pdbout,"TER\n");
}

typedef struct {
	char	*res;
	float	*acc, *bcc;
	Vec	*ca, *cb, *cc, *cd;
	int	len;
} Seq;

float	recp_area();

#define N 2

int	seqmat[NACID][NACID];
char	pdbcode1[55];

typedef struct {
	int	id, a, b;
	float	w;
} Sheet;

netsort (ac,bc) const void *ac, *bc;
{
Sheet   *a = (Sheet*)ac, *b = (Sheet*)bc;
        if (a->w < b->w) return 1;
        if (a->w > b->w) return -1;
        return 0;
}

pairsort (ac,bc) const void *ac, *bc;
{
Pairs   *a = (Pairs*)ac, *b = (Pairs*)bc;
        if (a->c < b->c) return 1;
        if (a->c > b->c) return -1;
        return 0;
}

lineup (code,prot,chain,secstr,gapen,betaw,betae,betaij)
char	*code;
Pdbentry_ *prot;
int	chain;
Secs	*secstr;
float	gapen, betaw, betae;
Pairs	*betaij;
{
Seq	seq[1];
int	n, i, j, k, len, lens;
float	**nets, *dens, *vars;
int	*ends;
FILE	*out, *ras, *var, *rod;
float	s = 0.15;
	k = 0;
	Ps(prot->Compound) NL
	len = protin(prot,chain,seq+k);
	for (i=0; i<=len; i++) vcopy(seq[k].ca[i], seq[k].cc+i);
	nets = (float**)alloca(sizeof(float*)*(len+4));
	for (i=0; i<len+4; i++) {
		nets[i] = (float*)alloca(sizeof(float)*(len+4));
		for (j=0; j<len+4; j++) nets[i][j] = 0.0;
	}
	ends = (int*)alloca(sizeof(int)*(len+1));
	dens = (float*)alloca(sizeof(float)*(len+1));
	vars = (float*)alloca(sizeof(float)*(len+1));
/* run here with a=1.65 b=1.2 e=0.23 */
	smooth(seq[k].ca,len,2);
	beta(seq[k].ca,len,ends,nets,betae);
/* run here with a=1.65 b=1.5 e=0.33 smooth(seq[k].ca,len,2); */
	out = fopen("test.out","w");
	rod = fopen("test.pdb","w");
	for (i=1; i<len+1; i++) { float d = nets[i][i+1] + nets[i][i-1];
                fprintf(out,"ATOM%7d  CA  GLY A%4d     %7.3f %7.3f %7.3f   0.00 %5.2f\n",
                        i, i, s*seq[k].ca[i].x, s*seq[k].ca[i].y, s*seq[k].ca[i].z, d);
		dens[i] = d;
	}
	fprintf(out,"TER\n");
	lens = lines(gapen,betaw,seq[k].ca+1,seq[k].cc+1,seq[k].cb,len,ends,dens,vars);
	for (i=1; i<=lens*2; i++) { float d = dens[i] - 1.0;
		d = 1.0+2.0*exp(-d*d*0.5);
                fprintf(out,"ATOM%7d  CA  GLY B%4d     %7.3f %7.3f %7.3f %5.2f %5.2f\n",
                        i, i, s*seq[k].cb[i].x, s*seq[k].cb[i].y, s*seq[k].cb[i].z, dens[i], d);
	}
	fprintf(out,"TER\n");
	for (i=1; i<=lens*2; i++) { float d = dens[i] - 1.0;
		d = 1.0+2.0*exp(-d*d*0.5);
                fprintf(rod,"ATOM%7d  CA  GLY B%4d     %7.3f %7.3f %7.3f %5.2f %5.2f\n",
                        i, i, seq[k].cb[i].x, seq[k].cb[i].y, seq[k].cb[i].z, dens[i], d);
	}
	fprintf(rod,"TER\n");
	var = fopen("vars.dat","w");
	ras = fopen("test.ras","w");
	fprintf(ras,"wireframe off\n");
	fprintf(ras,"select *A\n");
	fprintf(ras,"backbone 15\n");
	fprintf(ras,"colour grey\n");
	for (i=1; i<=lens; i++) { float pitch = dens[i*2], d = 0.0; char t='-';
		for (j=ends[i*2-1]; j<=ends[i*2]; j++)
			d += nets[j][j-1] + nets[j][j] + nets[j][j+1];
		fprintf(ras,"select %d-%dB\n", i*2-1,i*2);
		if (d>0.001 && pitch>2.0 ) {
			t = 'B'; fprintf(ras,"backbone  60\n");
		} else {
			if (pitch>2.20) {
				t = 'E'; fprintf(ras,"backbone  30\n");
			} else {
				if (pitch<1.75) { t = 'A'; fprintf(ras,"backbone 120\n"); }
					else       { t = '3'; fprintf(ras,"backbone  90\n"); }
			}
		}
		fprintf(ras,"colour temperature\n");
		for (j=ends[i*2-1]; j<=ends[i*2]; j++) secstr[j].seg = t;
		secstr[ends[i*2-1]].seg = tolower(t);
		secstr[ends[i*2]].seg = tolower(t);
		if (strlen(code)==4) { code[4]=' '; code[5] = (char)0; }
		if (vars[i*2]<0.0) vars[i*2] = 0.0;
		fprintf(var,"%c %s%5d%5d %7.4f %7.5f %3d\n",
			t, code, ends[i*2-1], ends[i*2], dens[i*2], sqrt(vars[i*2]), ends[i*2]-ends[i*2-1]+1);
	}
	n = 0;
	for (i=0; i<=len; i++) {
		for (j=0; j<=len; j++) {
			if (i<j) continue;
			if (i-j < 3) continue;
			if (nets[i][j] < NOISE) continue;
			Pi(i) Pi(j) Pr(nets[i][j]) NL
			betaij[n].a = i;
			betaij[n].b = j;
			betaij[n].c = nets[i][j];
			n++;
		}
	}
	betaij[n].c = -1.0;
	qsort(betaij,n,sizeof(Pairs),pairsort);
}

protin (prot,chain,seq)
Pdbentry_ *prot;
int	chain;
Seq *seq;
{
FILE	*pdb;
int	i, j, len;
	len = copyca(prot->Chains+chain,seq,0,1.0);
	add_cb(seq);
	return len;
}

add_cb (seq)
Seq	*seq;
{       int     i;
        for (i=1; i<=seq->len; i++)
        { Vec   n, c, b;
          float d, bond = 3.0;
                vsub(seq->ca[i],seq->ca[i-1],&n);
                vnorm(&n);
                vsub(seq->ca[i],seq->ca[i+1],&c);
                vnorm(&c);
                vadd(n,c,&b);
                vnorm(&b);
                vmul(&b,bond);
                vadd(seq->ca[i],b,&seq->cb[i]);
        }
}

extend (res,i,j,k,new)
Vec	*res;
int	i, j, k, new;
{
	Vec	m, v;
	vave(res[j],res[k],&m);
	vsub(m,res[i],&v);
	vadd(m,v,&res[new]);
}
 
copyca (pdb,s,flip,z)
Chain_  *pdb;
Seq	*s;
int	flip;
float	z;
{	int	i, n;
	char	*seq;
	Vec	*ca, *cb, *cc, *cd;
	float	*acc, *bcc;
	n = pdb->Aano;
	seq = (char*)malloc(sizeof(char)*(n+3));
	acc = (float*)malloc(sizeof(float)*(n+3));
	bcc = (float*)malloc(sizeof(float)*(n+3));
	ca = (Vec*)malloc(sizeof(Vec)*(n+3));
	cb = (Vec*)malloc(sizeof(Vec)*(n+3));
	cc = (Vec*)malloc(sizeof(Vec)*(n+3));
	cd = (Vec*)malloc(sizeof(Vec)*(n+3));
	for (i=0; i<n; i++) {
		ca[i+1].x = pdb->Atoms[i].X;
		ca[i+1].y = pdb->Atoms[i].Y;
		ca[i+1].z = pdb->Atoms[i].Z;
		acc[i+1] = pdb->Atoms[i].Occu;
		bcc[i+1] = pdb->Atoms[i].Bfact;
		seq[i+1] = pdb->Atoms[i].Aa;
		if (seq[i+1]<'A' || seq[i+1]>'Z') {
			printf("*NB* funny aa = %c\n", seq[i+1]);
			seq[i+1] = 'X';
		}
	}
	seq[0] = 'n';
        extend(ca,3,2,1,0);    
        extend(ca,n-2,n-1,n,n+1);
	seq[n+1] = 'c';
	seq[n+2] = 0;
	for (i=0; i<=n+1; i++) ca[i].z *= z;
	if (flip) flipseq(ca,seq,acc,n);
	s->res = seq;
	s->acc = acc;
	s->bcc = bcc;
	s->ca = ca;
	s->cb = cb;
	s->cc = cc;
	s->cd = cd;
	s->len = n;
	return n;
}

flipseq (ca,seq,acc,n) Vec *ca; char *seq; float *acc; int n;
{
int	i;
	for (i=0; i<=n/2; i++)
	{ Vec r; char c; float a;
	  int j = n+1-i;
		r = ca[i]; ca[i] = ca[j]; ca[j] = r;
		c = seq[i]; seq[i] = seq[j]; seq[j] = c;
		a = acc[i]; acc[i] = acc[j]; acc[j] = a;
	}
}
 
getca (res,pdb)
Vec    *res;
FILE	*pdb;
{	int	i = 1;
	char	line[225], junk[30];
        while(!feof(pdb)) {
		read_line(pdb,line);
		if (!strncmp(line,"TER",3)) break;
		if (strncmp(line,"ATOM",4)) continue;
		if (strncmp(line+13,"CA ",3)) continue;
		sscanf(line,"%30c %f%f%f",
                       	junk, &res[i].x, &res[i].y, &res[i].z);
		i++;
	}
	i--;
        extend(res,3,2,1,0);    
        extend(res,i-2,i-1,i,i+1);
	return i;
}
 
putpdb (seq,out,id)
Seq	*seq;
FILE    *out;
char	id;
{       int     i = 0, n = 0;
	int     len = seq->len;
        for (i=1; i<=len; i++) {
                fprintf(out,"ATOM%7d  CA  GLY%c%5d     %7.3f %7.3f %7.3f  0.00  0.00\n",
                        i, id, i, seq->ca[i].x, seq->ca[i].y, seq->ca[i].z);
        }
        fprintf(out,"TER\n");
}

setframe (a, b, c, frame)
    Vec a, b, c;
    Mat *frame;
{
    int    i;
    Vec    x, y, z ;
	vsub(c,a,&x);
	vave(c,a,&c);
	vsub(c,b,&y);
	vprod(y,x,&z);
	vprod(z,x,&y);
	vnorm(&x);
	vnorm(&y);
	vnorm(&z);
	VtoM(x,y,z,frame);
}

norm (sigcut, data, n) float sigcut, *data; int n; {
float   d, fn, ave, var, sig;
int     i, j, k, mods;
        fn = (float)n;
        ave = 0.0;
        for (i=0; i<n; i++) ave += data[i];
        ave /= fn;
        var = 0.0;
        for (i=0; i<n; i++) {
                var += data[i]*data[i];
        }
        var /= fn;
        sig = sqrt(var);
        mods = 0;
        for (i=0; i<n; i++) {
                data[i] /= sig;
                if (data[i] > sigcut) {
                        data[i] =  sigcut + 0.5*(data[i]-sigcut);
                        mods++;
                }
        }
        return mods;
}

normn (sigcut, data, m, n) float sigcut, **data; int m, n; {
float   d, fn, dmax,
        ave, var, sig;
int     i, j, k, mods;
        fn = (float)(m*n);
        ave = 0.0;
        for (i=1; i<=n; i++) for (j=1; j<=m; j++) ave += data[i][j];
        ave /= fn;
        var = 0.0;
        for (i=1; i<=n; i++) for (j=1; j<=m; j++) var += data[i][j]*data[i][j];
        var /= fn;
        sig = sqrt(var);
        mods = 0;
	dmax = 1.0;
        for (i=1; i<=n; i++) for (j=1; j<=m; j++) {
                data[i][j] /= sig;
                if (data[i][j] > sigcut) {
                        data[i][j] =  sigcut + 0.5*(data[i][j]-sigcut);
                        mods++;
                }
		if (data[i][j] > dmax) dmax = data[i][j];
        }
        for (i=1; i<=n; i++) for (j=1; j<=m; j++) data[i][j] /= dmax;
        return mods;
}

matin()
{
/*
Mutation Data Matrix (120 PAMs) + 8
*/
int	mat[529] = {
 3,-3,-1, 0,-3,-1, 0, 1,-3,-1,-3,-2,-2,-4, 1, 1, 1,-7,-4, 0, 0,-1,-1,
-3, 6,-1,-3,-4, 1,-3,-4, 1,-2,-4, 2,-1,-5,-1,-1,-2, 1,-5,-3,-2,-1,-2,
-1,-1, 4, 2,-5, 0, 1, 0, 2,-2,-4, 1,-3,-4,-2, 1, 0,-4,-2,-3, 3, 0,-1,
 0,-3, 2, 5,-7, 1, 3, 0, 0,-3,-5,-1,-4,-7,-3, 0,-1,-8,-5,-3, 4, 3,-2,
-3,-4,-5,-7, 9,-7,-7,-4,-4,-3,-7,-7,-6,-6,-4, 0,-3,-8,-1,-3,-6,-7,-4,
-1, 1, 0, 1,-7, 6, 2,-3, 3,-3,-2, 0,-1,-6, 0,-2,-2,-6,-5,-3, 0, 4,-1,
 0,-3, 1, 3,-7, 2, 5,-1,-1,-3,-4,-1,-3,-7,-2,-1,-2,-8,-5,-3, 3, 4,-1,
 1,-4, 0, 0,-4,-3,-1, 5,-4,-4,-5,-3,-4,-5,-2, 1,-1,-8,-6,-2, 0,-2,-2,
-3, 1, 2, 0,-4, 3,-1,-4, 7,-4,-3,-2,-4,-3,-1,-2,-3,-3,-1,-3, 1, 1,-2,
-1,-2,-2,-3,-3,-3,-3,-4,-4, 6, 1,-3, 1, 0,-3,-2, 0,-6,-2, 3,-3,-3,-1,
-3,-4,-4,-5,-7,-2,-4,-5,-3, 1, 5,-4, 3, 0,-3,-4,-3,-3,-2, 1,-4,-3,-2,
-2, 2, 1,-1,-7, 0,-1,-3,-2,-3,-4, 5, 0,-7,-2,-1,-1,-5,-5,-4, 0,-1,-2,
-2,-1,-3,-4,-6,-1,-3,-4,-4, 1, 3, 0, 8,-1,-3,-2,-1,-6,-4, 1,-4,-2,-2,
-4,-5,-4,-7,-6,-6,-7,-5,-3, 0, 0,-7,-1, 8,-5,-3,-4,-1, 4,-3,-5,-6,-3,
 1,-1,-2,-3,-4, 0,-2,-2,-1,-3,-3,-2,-3,-5, 6, 1,-1,-7,-6,-2,-2,-1,-2,
 1,-1, 1, 0, 0,-2,-1, 1,-2,-2,-4,-1,-2,-3, 1, 3, 2,-2,-3,-2, 0,-1,-1,
 1,-2, 0,-1,-3,-2,-2,-1,-3, 0,-3,-1,-1,-4,-1, 2, 4,-6,-3, 0, 0,-2,-1,
-7, 1,-4,-8,-8,-6,-8,-8,-3,-6,-3,-5,-6,-1,-7,-2,-6,12,-2,-8,-6,-7,-5,
-4,-5,-2,-5,-1,-5,-5,-6,-1,-2,-2,-5,-4, 4,-6,-3,-3,-2, 8,-3,-3,-5,-3,
 0,-3,-3,-3,-3,-3,-3,-2,-3, 3, 1,-4, 1,-3,-2,-2, 0,-8,-3, 5,-3,-3,-1,
 0,-2, 3, 4,-6, 0, 3, 0, 1,-3,-4, 0,-4,-5,-2, 0, 0,-6,-3,-3, 4, 2,-1,
-1,-1, 0, 3,-7, 4, 4,-2, 1,-3,-3,-1,-2,-6,-1,-1,-2,-7,-5,-3, 2, 4,-1,
-1,-2,-1,-2,-4,-1,-1,-2,-2,-1,-2,-2,-2,-3,-2,-1,-1,-5,-3,-1,-1,-1,-2
};
        int     n, i, j, m = 8;
        char    acid[24], c;
	strcpy(acid,"ARNDCQEGHILKMFPSTWYVBZX");
        printf("%s\n",acid);
        printf("matrix constant = %d\n",m);
	for (i=0; i<NACID; i++) for (j=0; j<NACID; j++) seqmat[i][j] = m;
	n = 0;
        for (i = 0; acid[i]; i++ )
        {       int ai = acid[i]-'A';
        	for (j = 0; acid[j]; j++ )
                {       int aj = acid[j]-'A';
                        seqmat[ai][aj] = mat[n] + m;
			n++;
                }
        }
}

oldmatin(file,mat)
        char    *file;
        int     mat[NACID][NACID];
{
        int     i, j, mat_const;
        char    acid[NACID], c;
        FILE    *mat_file;

        mat_file = fopen(file,"r");
        while( c = getc(mat_file), c != '\n' ) putchar(c); NL
        fscanf(mat_file,"%s\n",acid);
        printf("%s\n",acid);
        fscanf(mat_file,"%d\n",&mat_const);
        printf("matrix constant = %d\n",mat_const);
        for( i = 0; acid[i]; i++ ) 
        {       int     ai = acid[i]-'A';
                for( j = 0; acid[j]; j++ ) 
                {       int aj = acid[j]-'A';
                        fscanf(mat_file,"%d",&mat[ai][aj]);
                        mat[ai][aj] += mat_const;
                }
        }
}

#define TURN 3
#define BEND 6

beta (ca,n,ends,net,betae) Vec *ca; int n, *ends; float **net, betae;
{
int	i, j, k, m, nn, min, max;
int	*mates[2];
Sheet	*sheet;
float   **mat, *ave, up = betae;
float	cu = 4.7+up, cut = cu+up, cutt = cut+up;
	ave = (float*)alloca(sizeof(float)*(n+2));
	mat = (float**)alloca(sizeof(float*)*(n+2));
	mates[0] = (int*)alloca(sizeof(int)*(n+2));
	mates[1] = (int*)alloca(sizeof(int)*(n+2));
	sheet = (Sheet*)alloca(sizeof(Sheet)*(n+2));
	for (i=0; i<n+2; i++) {
	        mat[i] = (float*)alloca(sizeof(float)*(n+2));
		for (j=0; j<n+2; j++) net[i][j] = 0.0;
		mates[0][i] = mates[1][i] = -1;
		sheet[i].w = -1.0;
	}
	for (i=1; i<n; i++) for (j=i+1; j<=n; j++) mat[i][j] = mat[j][i] = vdif(ca[i],ca[j]);
	for (i=2; i<n; i++)
	{ int	p[3], q[3], r, s[3];
	  float dmin;
		p[0]=p[1]=p[2]=q[0]=q[1]=q[2] = 0;
		s[0]=i-1, s[1]=i, s[2]=i+1;
		dmin = 999.9;
		for (j=2; j<n; j++)
		{ float	a = mat[i][j];
			if (abs(i-j)<TURN) continue;
			if (a > cut) continue;
			if (a > dmin) continue;
			p[1]=j; dmin=a;
		}
		if (p[1] && mat[s[1]][p[1]] < cu ) {
			if (mat[s[0]][p[1]-1]<cut) p[0] = p[1]-1;
			if (mat[s[0]][p[1]+1]<cut) p[0] = p[1]+1;
			if (mat[s[2]][p[1]-1]<cut) p[2] = p[1]-1;
			if (mat[s[2]][p[1]+1]<cut) p[2] = p[1]+1;
			if (p[0] && p[2] && abs(s[2]-p[2])>TURN
					 && abs(s[0]-p[0])>TURN ) {
				net[s[1]][s[0]] += 1.0;
				net[s[0]][s[1]] += 1.0;
				net[s[1]][s[2]] += 1.0;
				net[s[2]][s[1]] += 1.0;
				net[s[1]][p[1]] += 1.0;
				net[p[1]][s[1]] += 1.0;
			}
		}
		dmin = 999.9;
		for (j=2; j<n; j++)
		{ float	a = mat[i][j];
			if (abs(i-j)<TURN) continue;
			if (a > cut) continue;
			if (a > dmin) continue;
			if (abs(j-p[1])<BEND) continue;
			q[1]=j; dmin=a;
		}
		if (!p[1] || !q[1]) continue;
		if ((mat[s[0]][p[1]-1] < mat[s[0]][p[1]+1])
		 && (mat[s[2]][p[1]+1] < mat[s[2]][p[1]-1])) {
			p[0]=p[1]-1; p[2]=p[1]+1;
		} else {
			p[0]=p[1]+1; p[2]=p[1]-1;
		}
		if ((mat[s[0]][q[1]-1] < mat[s[0]][q[1]+1])
		 && (mat[s[2]][q[1]+1] < mat[s[2]][q[1]-1])) {
			q[0]=q[1]-1; q[2]=q[1]+1;
		} else {
			q[0]=q[1]+1; q[2]=q[1]-1;
		}
		if (mat[s[0]][p[0]] > cutt) continue;
		if (mat[s[2]][p[2]] > cutt) continue;
		if (mat[s[0]][q[0]] > cutt) continue;
		if (mat[s[2]][q[2]] > cutt) continue;
		if (abs(s[2]-p[2])<TURN) continue;
		if (abs(s[0]-p[0])<TURN) continue;
		if (abs(s[2]-p[0])<TURN) continue;
		if (abs(s[0]-p[2])<TURN) continue;
		if (abs(s[2]-q[2])<TURN) continue;
		if (abs(s[0]-q[0])<TURN) continue;
		if (abs(s[2]-q[0])<TURN) continue;
		if (abs(s[0]-q[2])<TURN) continue;
		if (abs(q[2]-p[2])<TURN) continue;
		if (abs(q[0]-p[0])<TURN) continue;
		if (abs(q[2]-p[0])<TURN) continue;
		if (abs(q[0]-p[2])<TURN) continue;
		printf("sheet: %4d %4d %4d\n", p[1],s[1],q[1]);
		mates[0][s[0]] = p[0];
		mates[0][s[1]] = p[1];
		mates[0][s[2]] = p[2];
		mates[1][s[0]] = q[0];
		mates[1][s[1]] = q[1];
		mates[1][s[2]] = q[2];
		/* link across chain */
		net[s[0]][q[0]] += 1.0; net[q[0]][s[0]] += 1.0;
		net[s[0]][p[0]] += 1.0; net[p[0]][s[0]] += 1.0;
		net[s[1]][q[1]] += 1.0; net[q[1]][s[1]] += 1.0;
		net[s[1]][p[1]] += 1.0; net[p[1]][s[1]] += 1.0;
		net[s[2]][q[2]] += 1.0; net[q[2]][s[2]] += 1.0;
		net[s[2]][p[2]] += 1.0; net[p[2]][s[2]] += 1.0;
		/* link along chain */
		net[s[1]][s[0]] += 1.0; net[s[0]][s[1]] += 1.0;
		net[s[1]][s[2]] += 1.0; net[s[2]][s[1]] += 1.0;
		net[p[1]][p[0]] += 1.0; net[p[0]][p[1]] += 1.0;
		net[p[1]][p[2]] += 1.0; net[p[2]][p[1]] += 1.0;
		net[q[1]][q[0]] += 1.0; net[q[0]][q[1]] += 1.0;
		net[q[1]][q[2]] += 1.0; net[q[2]][q[1]] += 1.0;
		/* extra links along chain 
		net[s[0]][s[0]-1] += 1.0; net[s[0]-1][s[0]] += 1.0;
		net[s[2]][s[2]+1] += 1.0; net[s[2]+1][s[2]] += 1.0;
		*/
	}
	NL
	nn = 0;
	for (i=1; i<=n; i++)
	{ float max1, max2, sum;
	  int	mate1, mate2;
		mate1 = mate2 = -1;
		max1 = max2 = sum = 0.0;
		for (j=1; j<=n; j++) { float nij;
			if (abs(i-j) < 2) continue;
			nij = net[i][j];
			if (nij<0.5) continue;
			sum += nij;
			if (nij > max2) {
				if (nij > max1) {
					mate2 = mate1; mate1 = j;
					max2 = max1; max1 = nij;
				} else {
					mate2 = j; max2 = nij;
				}
				if (abs(mate1-mate2) < 4) {
					mate2 = -1; max2 = 0.0;
				}
			}
		}
		if (mate1>0 && mate2>0 && (mate1*mate2 != mates[0][i]*mates[1][i])) {
			Pi(i) Pi(mate1) Pi(mate2) Pi(mates[0][i]) Pi(mates[1][i]) NL
			if ((mates[0][i]<0) && (mates[1][i]<0)) continue;
		}
		if (sum<0.5) continue;
		sheet[nn].id = i;
		sheet[nn].w = sum;
		sheet[nn].a = mate1;
		sheet[nn].b = mate2;
		nn++;
	}
	qsort(sheet,nn,sizeof(Sheet),netsort);
	for (i=0; i<nn; i++)
	{ int id, a, b;
	  Vec s, t;
	  float f;
		a = sheet[i].a;
		b = sheet[i].b;
		if (b<0 || a<0) continue;
		id = sheet[i].id;
		vave(ca[a],ca[b], &t);
		vsub(t,ca[id], &s);
		f = sheet[i].w/12.0;
		vmul(&s,f);
		vadd(s,ca[id], ca+id);
	}
	for (i=0; i<=n; i++) ends[i] = 0;
	for (i=3; i<n-2; i++)
	{ float lasts = net[i-3][i-2]+net[i-2][i-1]+net[i-1][i], next = net[i][i+1],
		nexts = net[i+3][i+2]+net[i+2][i+1]+net[i+1][i], last = net[i][i-1];
		if ((lasts<0.5) && (next>0.5)) ends[i-1] = 1;
		if ((nexts<0.5) && (last>0.5)) ends[i+1] = 1;
	}
	min = 9999; max = 0;
	for (i=2; i<n; i++) {
		if (net[i][i+1]>0.0001 && i<min) min = i-1;
		if (net[i-1][i]>0.0001 && i>max) max = i+1;
	}
	if (min<9999) ends[min-1] = 1;
	if (max>0 )   ends[max+1] = 1;
}

consec_seen (chain,secs,pcta,pctb,pctc,secstr)
Chain_ *chain;
char *secs;
float	*pcta, *pctb, *pctc;
Secs	*secstr;
{
int	breaks = 0;
int	i, j, len = chain->Aano;
float alpha[9] = {0.0,  3.82,  5.41,  5.04,  6.21,  8.66,  9.82, 10.52, 12.37 };
float fre10[9] = {0.0,  3.80,  5.40,  5.90,  8.50, 10.50, 12.50, 14.00, 16.00 };
float betas[9] = {0.0,  3.75,  6.89, 10.44, 13.78, 17.29, 20.67, 24.16, 27.56 };
float noise = 0.35;
float     *psec, **mat;
	psec = (float*)alloca(sizeof(float)*(len+2));
	mat = (float**)alloca(sizeof(float*)*(len+2));
        for (i=1; i<=len; i++) {
		mat[i] = (float*)alloca(sizeof(float)*(len+2));
        	for (j=1; j<=len; j++) {
			mat[i][j] = cadist(chain->Atoms,i-1,j-1);
		}
		if (i==1) continue;
		if (mat[i][i-1] > 4.0) breaks++;
	}
	*pcta = *pctb = 0.0;
        for (i=1; i<=len; i++)
        { float rmsa, rmsb, rmsc;
                rmsa = 1.15 * dmatch(mat,i,alpha,3,len);
                rmsb = 1.35 * dmatch(mat,i,betas,3,len);
                rmsc = 0.85 * dmatch(mat,i,fre10,3,len);
                secs[i] = '-';
                psec[i] = -1.0;
                if (rmsa>noise && rmsa>rmsb && rmsa>rmsc) {
                        secs[i] = 'H';
                        psec[i] = rmsa;
                }
                if (rmsb>noise && rmsb>rmsa && rmsb>rmsc) {
                        secs[i] = 'E';
                        psec[i] = rmsb;
                }
                if (rmsc>noise && rmsc>rmsb && rmsc>rmsa) {
                        secs[i] = '3';
                        psec[i] = rmsc;
                }
		if (secs[i]=='H') { secstr[i].rms = 'A'; *pcta += 1.0; }
		if (secs[i]=='E') { secstr[i].rms = 'B'; *pctb += 1.0; }
		if (secs[i]=='3') { secstr[i].rms = '3'; *pctc += 1.0; }
        }
	*pcta = *pcta * 100.0/(float)len;
	*pctb = *pctb * 100.0/(float)len;
	*pctc = *pctc * 100.0/(float)len;
	return breaks;
}

float cadist (a,i,j)
Atom_ *a;
int	i, j;
{
float	x = a[i].X - a[j].X,
	y = a[i].Y - a[j].Y,
	z = a[i].Z - a[j].Z;
	return sqrt(x*x + y*y + z*z);
}
	
float dmatch (mat, m, str, n, l)
float     **mat;
int     m, n, l;
float   *str;
{
float   score, sum, d, in, nn = (float)(n+n+1);
int     i, j;
        nn = nn*nn-nn;
        in = sum = 0.0;
        for (i=-n; i<=n; i++) {
                for (j=-n; j<=n; j++)
                { int   mi = m+i, mj = m+j;
                        if (i==j) continue;
                        if (mi<1 || mj<1) continue;
                        if (mi>l || mj>l) continue;
                        d = mat[mi][mj] - str[abs(i-j)];
                        sum += d*d;
                        in += 1.0;
                }
        }
        score = nn/(nn+sum);
        score = score*in/nn;
        return score;
}

int axisfit(a,m,n,out) Vec *a; int m,n; Vec *out;
{
float	ratio[22] = { 0.1,1.0/10.0,2.0/10.0,3.0/10.0,4.0/10.0,5.0/10.0,6.0/10.0,7.0/10.0,8.0/10.0,9.0/10.0,  // flat
                          1.0,10.0/9.0,10.0/8.0,10.0/7.0,10.0/6.0,10.0/5.0,10.0/4.0,10.0/3.0,10.0/2.0,10.0}; // long
int     i, j, in = 0;
float   rog, **mom, **vec, val[3];
Vec   	c, cog, axis;
float	r, ab, bc, fit;
int	rat;
	mom = (float**)alloca(sizeof(float*)*3);
	vec = (float**)alloca(sizeof(float*)*3);
        for (i=0; i<3; i++) {
                mom[i] = (float*)alloca(sizeof(float)*3);
                vec[i] = (float*)alloca(sizeof(float)*3);
        }
        vinit(&cog);
        for (i=m; i<=n; i++) { vsum(a[i],&cog); in++; }
        vdiv(&cog,(float)in);
	Pi(m) Pi(n) Pi(in) Pv(cog) NL
	vcopy(cog,out); // return centre of gravity
        for (i=0; i<3; i++) for (j=0; j<3; j++) mom[i][j] = 0.0;
	rog = 0.0;
        for (i=m; i<=n; i++) {
                vsub(a[i],cog,&c);
		rog += vsqr(c);
                mom[0][0] += c.x * c.x;
                mom[1][1] += c.y * c.y;
                mom[2][2] += c.z * c.z;
                mom[0][1] = mom[1][0] += c.x * c.y;
                mom[0][2] = mom[2][0] += c.x * c.z;
                mom[1][2] = mom[2][1] += c.y * c.z;
        }
	rog /= (float)in;
	Pr(rog) NL
	project(mom,3,val,vec);
	printf(" Eigen values = %f %f %f\n",val[0],val[1],val[2]); 
	for (i=0; i<3; i++) val[i] = sqrt(val[i]);
	// return sorted unit axes and length (of equivalent ellipsoid)
	axis.x = vec[0][0]; axis.y = vec[0][1]; axis.z = vec[0][2]; vcopy(axis,out+1);
	axis.x = vec[1][0]; axis.y = vec[1][1]; axis.z = vec[1][2]; vcopy(axis,out+2);
	axis.x = vec[2][0]; axis.y = vec[2][1]; axis.z = vec[2][2]; vcopy(axis,out+3);
	axis.x = val[0]; axis.y = val[1]; axis.z = val[2]; vcopy(axis,out+4);
	ab = axis.y/axis.x;
	bc = axis.z/axis.y;
	if (bc > ab) { // smallest pair are close so prolate (cigar) 
		r = axis.x/((axis.y+axis.z)*0.5);
	} else {	// biggest pair are closer so oblate (flying saucer)
		r = axis.z/((axis.y+axis.x)*0.5);
	}
	fit = 999999.9;
	for (i=1; i<20; i++) { float d;
		d = ratio[i]-r; d*=d;
		if (d < fit) { fit = d; rat = i; }
	}
	return rat;
}
