/*
cc -O sap.c -o sap util/aa/cones.o util/aa/stutest.o util/aa/bestrot.o util/wt/util.o util/wt/geom.o util/wt/sort.o util/aa/pdbprot.o util/aa/matrix.o util/aa/siva.o util/aa/ql.o util/eigen.o -lm

*/
#include <alloca.h>
#include "util/wt/incl/util.h"
#include "util/wt/incl/geom.h"
#include "util/aa/incl/pdbprot.h"
#include "util/aa/incl/matrix.h"
#define NALLOC 1000
#define NACID 30

typedef struct {
	Vec	v;
	float	d;
	Vec	cos;
} Tri;

typedef struct {
	int	c, r;
} Cel;

typedef struct {
	int	a, b;
	float	c;
} Pairs;

typedef struct {
	char	*res;
	float	*acc;
	Vec	*ca, *cb;
	int	len;
} Seq;

float	compare();
float	recycle();
float	score();
float	rescore();
float	match_list();
float	add_path();
float	get_path();
float	local_rms();
float	super();
double	supermac();
double  drand48();

#define CYCLES 5
#define BIAS_WT 0.05
#define BIAS_DAMP 0.5
#define PATHSUM_WT 0.0
#define PATH_WT 1.0
#define N 2

int	self;
int	seqmat[NACID][NACID];

char	pdbcode[99];

main(argc,argv)
int argc; char *argv[];
{
Tri	**a[N];
Cel	**cel[N];
Seq	seq[N];
char	id[N*4+1], file1[255], file2[255];
float	**scores;
int	i, j, cycles;
char	*c;
float	rms, z = 1.0;
Pdbentry_ *prot1, *prot2;
long    rseed = 234;
        srand48(rseed);
	if (argc < 3) {
		printf("usage:  sap file1.pdb file2.pdb\n");
		exit(1);
	} else {
		strcpy(file1,argv[1]);
		strcpy(file2,argv[2]);
	}
	strcpy(pdbcode,argv[2]);
	Ps(file1) NL
	prot1 = get_pdb(file1,1,1);
	Ps(prot1->Compound) NL
	cones(prot1);
	Ps(file2) NL
	prot2 = get_pdb(file2,1,1);
	Ps(prot2->Compound) NL
	cones(prot2);
	matin();
	protin(prot1,seq+0,a+0,cel+0,1.0,0);
	protin(prot2,seq+1,a+1,cel+1,1.0,0);
	self = 0;
	if (!strcmp(file1,file2)) self = 1;
	cycles = CYCLES;
	if (argc==4) {
		sscanf(argv[3],"%d", &cycles);
                printf("\n%d cycles\n", cycles);
	}
	rms = compare(cycles,a[1],a[0],cel[1],cel[0],seq+1,seq+0,1);
}

stats (half, data, n) float **data; int half, n; {
double	fn1, fn2, ave1, ave2, var1, var2;
float	sig1, sig2, dmax, dmin, smax, clear, score, noise;
int     i, j, k, n1, n2;
	n1 = n2 = 0;
        ave1 = ave2 = 0.0;
	smax = dmax = 0.0;
        for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			if (half && i<j) continue;
			k = i+j+1; 
			if (k-(2*(k/2))) continue;
			k = i+j-1;
			if ((i-j)*(i-j)==1 && !(k-(4*(k/4))) ) {
				if (dmax<data[i][j]) dmax = data[i][j];
				ave1 += data[i][j];
				n1++;
			} else {
				if (smax<data[i][j]) smax = data[i][j];
				ave2 += data[i][j];
				n2++;
			}
		}
	}
	fn1 = (double)n1; fn2 = (double)n2;
        ave1 /= fn1;
        ave2 /= fn2;
        var1 = var2 = 0.0;
	dmin = dmax;
	noise = 0.0;
        for (i=0; i<n; i++) { float d;
		for (j=0; j<n; j++) {
			if (half && i<j) continue;
			k = i+j+1;
			if (k-(2*(k/2))) continue;
			k = i+j-1;
			if ((i-j)*(i-j)==1 && !(k-(4*(k/4))) ) {
				if (dmin>data[i][j]) dmin = data[i][j];
				d = data[i][j] - ave1;
				var1 += d*d;
			} else {
				d = data[i][j] - ave2;
				var2 += d*d;
			}
			if (half || i<j) continue;
			d = data[i][j] - data[j][i];
			noise += d*d;
		}
	}
        var1 /= fn1; var2 /= fn2;
	stutest(ave1,ave2,var1,var2,n1,n2);
	sig1 = sqrt(var1); sig2 = sqrt(var2);
	noise = sqrt(2.0*noise/(fn1+fn2));
	score = (data[1][0] - ave2)/sig2;
	clear = (dmin-smax)/sig2;
	dmax = (dmax - ave2)/sig2;
	dmin = (dmin - ave2)/sig2;
	NLL
	printf("Alignment score = %5.3f StD(2) above mean controls\n", score);
	printf("Worst score = %5.3f StD(2) above best control\n", clear);
	printf("Max score = %5.3f, Min score = %5.3f\n", dmax, dmin);
	printf("StD(1) = %5.3f, StD(2) = %5.3f\n", sig1,sig2);
	if (half) return;
	printf("RMS alignment order noise = %5.3f\n\n", noise);
}

selsort (ac,bc) const void *ac, *bc;
{
Pairs   *a = (Pairs*)ac, *b = (Pairs*)bc;
        if (a->c < b->c) return 1;
        if (a->c > b->c) return -1;
        return 0;
}

float compare (cycles,a,b,c,d,seqa,seqb,print)
Tri	**a, **b; 
Cel	**c, **d; 
Seq	*seqa, *seqb;
int	cycles, print;
{
float	**bias, **sim, **sec, score;
int	i, j, m, n, cycle, nsel;
int	lena, lenb;
Pairs	*sel;
	lena = seqa->len;
	lenb = seqb->len;
Pi(lena) Pi(lenb) NL
        sel = (Pairs*)alloca(sizeof(Pairs)*lena*lenb);
        sec = (float**)alloca(sizeof(float*)*(lenb+2));
        sim = (float**)alloca(sizeof(float*)*(lenb+2));
        bias = (float**)alloca(sizeof(float*)*(lenb+2));
        for (i=0; i<lenb+2; i++) {
                sim[i] = (float*)alloca(sizeof(float)*(lena+2));
                sec[i] = (float*)alloca(sizeof(float)*(lena+2));
                bias[i] = (float*)alloca(sizeof(float)*(lena+2));
                for (j=0; j<lena+2; j++) bias[i][j] = 0.0;
        }
        for (i=3; i<lenb-1; i++) {
                for (j=3; j<lena-1; j++) {
		  float acc, rms, seq, aa, ab, dab;
		  int	ra, rb;
			if (i==j) sim[i][j] = bias[i][j] = 100.0; else continue;
			aa = seqa->acc[j];
			ab = seqb->acc[i];
			acc = aa-ab;
			acc = exp(-acc*acc);
			rms = local_rms(j,i,lena,lenb,a,b);
			rms = exp(-rms*rms);
			sec[i][j] = rms;
			ra = seqa->res[j] - 'A';
			rb = seqb->res[i] - 'A';
			seq = log(1.0+(float)seqmat[ra][rb]);
			bias[i][j] = acc+rms;
			if (print>1) { float d;
				d = vdif(seqa->ca[j],seqb->ca[i]);
				d = exp(-d*d*0.01)*(float)print;
				bias[i][j] += d;
			}
			bias[i][j] *= (float)cellhits(c[j],d[i],a,b,j,i,lena,lenb);
			sim[i][j] = bias[i][j];
                        if (self)
                        { float off;
                                if (i>=j) {
                                        bias[i][j] = 0.0;
                                } else {
                                        off = 4.0*(float)(j-i)/(float)lena - 2.0;
                                        bias[i][j] *= exp(-off*off);
                                }
                        }
		}
        }
	nsel = 0;
	recycle(0,seqa,seqb,bias,sel,sec,sim,a,b,&nsel,0);
/*
	for (cycle=1; cycle<cycles; cycle++) {
		if (print) printf("Cycle %d, %d residues selected\n", cycle, nsel);
		recycle(cycle,seqa,seqb,bias,sel,sec,sim,a,b,&nsel,0);
	}
*/
	cycle = cycles;
	if (print) printf("Cycle %d, %d residues selected\n", cycle, nsel);
	return  recycle(cycle,seqa,seqb,bias,sel,sec,sim,a,b,&nsel,print);
}

float local_rms(m,n,lena,lenb,a,b)
int     m,n, lena,lenb;
Tri     **a, **b;
{
float   sum, d;
int     i, j, k=3;
        sum = 0.0;
        for (i=-k; i<=k; i++) {
                if (m+i<=0 || n+i<=0) continue;
                if (m+i>=lena || n+i>=lenb) continue;
                for (j=-k; j<=k; j++) {
                        if (m+j<=0 || n+j<=0) continue;
                        if (m+j>=lena || n+j>=lenb) continue;
                        d = a[m+i][m+j].d - b[n+i][n+j].d;
                        sum += d*d;
                }
        }
        k = 2*k+1; k = k*k;
        return sum/(float)k;
}

float recycle (cycle,seqa,seqb,bias,sel,sec,sim,a,b,nsel,print)
int	cycle;
Seq	*seqa, *seqb;
float	**bias, **sec, **sim;
Tri	**a, **b; 
Pairs	*sel;
int	*nsel, print;
{
int	**aln, len, i, j, n;
int     lena, lenb, onaln;
float	score, acut, cyc_no = (float)cycle;
	lena = seqa->len;
	lenb = seqb->len;
	if (*nsel) score_pair(cyc_no,bias,sel,sim,a,b,lena,lenb,*nsel);
	aln = (int**)alloca(sizeof(int*)*2); TEST(aln)
	for (i=0; i < 2; i++) {
		aln[i] = (int*)alloca(sizeof(int)*(lena+lenb)); TEST(aln[i])
	}
	for (i=1; i<=lenb; i++) {
		for (j=1; j<=lena; j++) {
			if (i==j) sim[i][j] = 10.0; else continue;
			sim[i][j] += 0.1*bias[i][j];
                        if (self)
                        { float off;
                                if (i>=j) {
                                        sim[i][j] = 0.0;
                                } else {
                                        off = 4.0*(float)(j-i)/(float)lena - 2.0;
                                        sim[i][j] *= exp(-off*off);
                                }
                        }
		}
	}
/*
len = lenb;
for (i=1; i<=lenb; i++) aln[1][i] = aln[2][i] = i;
for (i=1; i<=lenb; i++) sim[i][i] = 100;
	score = get_path(aln,sim,lena,lenb,&len);
Pr(score) Pi(len) NL
*/
	if (print) {
		for (i=1; i<=lenb; i++) {
			aln[0][i] = aln[1][i] = i;
			sim[i][i] = 10.0;
		}
Pi(lenb) NL
		return super(a, b, seqa, seqb, sim, aln, lenb);
		exit(1);
	}
	if (print) {
		NL Pr(score) NLL
		onaln = check_sel(aln,sel,*nsel,len);
		printf("Percent sel on aln = %7.2f\n",
			100.0*(float)onaln/(float)*nsel);
		printf("Percent aln in sel = %7.2f\n",
			100.0*(float)onaln/(float)len);
		return super(a, b, seqa, seqb, sim, aln, len);
	}
/*
	for (i=1; i<=lenb; i++) {
		for (j=1; j<=lena; j++) {
			bias[i][j] *= BIAS_DAMP;
		}
	}
	for (i=len; i>0; i--) 
	{ int	p = aln[0][i],
		q = aln[1][i];
		bias[q][p] += log(1.0+sim[q][p]) * BIAS_WT;
	}
	normn(3.0,bias,lena,lenb);
        if (!cycle) { float sab, sa, sb;
                sa = (float)(aln[0][1]-aln[0][len]);
                sb = (float)(aln[1][1]-aln[1][len]);
                if (sa>sb) sab = sb/sa; else sab = sa/sb;
                sab *= 3.0;
                for (i=len; i>0; i--)
                { int   p=aln[0][i], q=aln[1][i];
                        bias[q][p] *= sab;
                }
        }
	n = 0;
	for (i=3; i<=lenb-1; i++) {
		for (j=3; j<=lena-1; j++) {
			sel[n].a = j;
			sel[n].b = i;
			sel[n].c = bias[i][j];
Pi(i) Pi(j) Pi(n) Pr(bias[i][j]) NL
			n++;
		}
        }
	if (n==0) { printf("*NB* No pairs selected\n"); exit(1); }
	qsort(sel,n-1,sizeof(Pairs),selsort);
	*nsel = 10+(int)(0.05*sqrt((float)(lena*lenb))*(1.0+cyc_no));
        if (self) *nsel /= 2;
	if (*nsel>100) *nsel = 100;
	return score;
*/
for (i=1; i<lenb; i++) sel[i].a = sel[i].b = i;
for (i=1; i<lenb; i++) sel[i].c = 100.0;
return 1.0;
}

cellhits (lista,listb,a,b,m,n,lena,lenb)
Cel	*lista, *listb;
Tri	**a, **b;
int	m, n, lena, lenb;
{
int	i, hits, na, nb, k, ka, kb;
	ka = kb = hits = 1;
	while (lista[ka].c>0 && listb[kb].c>0) { int ja, jb;
		if (lista[ka].c > listb[kb].c) { ka++; continue; }
		if (listb[kb].c > lista[ka].c) { kb++; continue; }
		/* same cell values */
		if (lista[ka].r==m) { ka++; continue; }
		if (listb[kb].r==n) { kb++; continue; }
		/* neither is the central res */
		for (ja=ka; lista[ja].c>0 && lista[ja].c==lista[ka].c; ja++) {
			for (jb=kb; listb[jb].c>0 && listb[jb].c==listb[kb].c; jb++) {
				if (lista[ja].r>m && listb[jb].r<n) continue;
				if (lista[ja].r<m && listb[jb].r>n) continue;
				hits++;
			}
		}
		ka = ja;
	}
	return hits;
}

check_sel (aln,sel,nsel,naln) Pairs *sel; int **aln, nsel, naln; {
int	i, j, n = 0;
	for (i=0; i<nsel; i++)
	{ int	sela = sel[i].a,
		selb = sel[i].b;
		for (j=1; j<=naln; j++)
		{ int	alna = aln[0][j],
			alnb = aln[1][j];
			if (sela==alna && selb==alnb) {
				n++;
				break;
			}
		}
	}
	return n;
}

trace_mat (mat,lena,lenb) int **mat, lena, lenb; {
int	i, j;
	NL Pi(lena) Pi(lenb) NL
	for (i=1; i<=lenb; i++) { 
		for (j=1; j<=lena; j++) printf("%3d", mat[i][j]); NL
	} NL
}

print_mat (scale, mat,lena,lenb) float scale, **mat; int lena, lenb; {
int	i, j;
	NL Pi(lena) Pi(lenb) NL
	for (i=1; i<=lenb; i++) { 
		for (j=1; j<=lena; j++)
		{ char	c;
	  	  float sij = mat[i][j]*scale;
			if (sij<10.0) c = '0'+(int)sij;
				else c = 'A'+(int)(sij*0.1)-1;
				if (sij > 260.0) c = '*';
				if (c=='0') c = '.';
			printf("%c", c);
		} NL
	} NL
}

score_pair (dif_wt,bias,sel,sim,a,b,la,lb,nsel)
float	dif_wt;
float	**bias, **sim;
Tri	**a, **b; 
Pairs	*sel;
int	la, lb, nsel;
{
int	i, j, k, l;
        for (i=0; i<lb+2; i++) {
                for (j=0; j<la+2; j++) sim[i][j] = 0.0;
        }
	for (i=0; i<nsel; i++) {
		score(dif_wt,bias,sim,a,b,sel[i].a,sel[i].b,la,lb);
	}
}

#define SEQW 0.001
#define CYCW 0.1
#define DISW 0.001
#define COSW 1.0
#define VECW 0.1

float score (wt,bias,sim,a,b,m,n,la,lb)
float	wt;
float	**bias, **sim;
Tri	**a, **b; 
int	m, n, la, lb;
{
int	i, j, k, na, nb, minlen;
float   **smn, path_score;
        smn = (float**)alloca(sizeof(float*)*(lb+2));
        for (i=0; i<lb+2; i++) {
                smn[i] = (float*)alloca(sizeof(float)*(la+2));
                for (j=0; j<la+2; j++) smn[i][j] = 0.0;
        }
	for (i=1; i<=lb; i++) {
                for (j=1; j<=la; j++)
		{ float cxd, cyd, czd, dif, gwt, swt, add, bdd, cdd, dwt, d;
			if (i<=n && j>=m) continue;
			if (i>=n && j<=m) continue;
                        d = (float)(abs(n-i)+abs(m-j));
                        swt = 1.0-exp(-d*d*SEQW);
			cxd = a[m][j].cos.x - b[n][i].cos.x; cxd = exp(-cxd*cxd*COSW);
			cyd = a[m][j].cos.y - b[n][i].cos.y; cyd = exp(-cyd*cyd*COSW);
			czd = a[m][j].cos.z - b[n][i].cos.z; czd = exp(-czd*czd*COSW);
			cdd = cxd*cyd*czd;
			dif = vddif(a[m][j].v,b[n][i].v);
			dif = exp(-dif*VECW);
			gwt = exp(-wt*wt*CYCW);
			d = a[m][j].d + b[n][i].d;
			dwt = exp(-d*d*DISW);
			smn[i][j] = 10.0*dwt*swt*dif*cdd + gwt*bias[i][j];
                        if (self && i>=j) smn[i][j] = 0.0;
		}
	}
	path_score = add_path(sim,smn,la,lb,m,n);
}

float rescore (a,b,m,n,aln,len)
Tri     **a, **b;
int     m, n, **aln, len;
{
        int     i, j, k;
        float   sum = 0.0;
                for (k=len; k>0; k--)
                { float cxd, cyd, czd, dif, gwt, swt, add, bdd, cdd, dwt,d;
                        i = aln[1][k];  j = aln[0][k];
                        d = (float)(abs(n-i)+abs(m-j));
                        swt = 1.0-exp(-d*d*SEQW);
/*
                        add = (float)(n-i);
                        bdd = (float)(m-j);
			d = add*add + bdd*bdd;
                        swt = 1.0-exp(-d*SEQW);
*/
                        cxd = a[m][j].cos.x - b[n][i].cos.x; cxd = exp(-cxd*cxd*COSW);
                        cyd = a[m][j].cos.y - b[n][i].cos.y; cyd = exp(-cyd*cyd*COSW);
                        czd = a[m][j].cos.z - b[n][i].cos.z; czd = exp(-czd*czd*COSW);
                        cdd = cxd*cyd*czd;
                        dif = vddif(a[m][j].v,b[n][i].v); dif = exp(-dif*VECW);
                        d = a[m][j].d + b[n][i].d;
                        dwt = exp(-d*d*DISW);
/*
                        add = a[m][j].d;
			bdd = b[n][i].d;
			d = add*add + bdd*bdd;
                        dwt = exp(-d*DISW);
*/
                        sum += 10.0*dwt*swt*dif*cdd;
                }
        return sum;
}

float add_path (sim,smn,na,nb,m,n) 
float	**sim, **smn;
int	na, nb, m, n;
{
int	**aln, len, i;
float	**s, score = 0.0;
	aln = (int**)alloca(sizeof(int*)*2); TEST(aln)
	for (i=0; i < 2; i++) {
		aln[i] = (int*)alloca(sizeof(int)*(na+nb)); TEST(aln[i])
	}
	if (m>1 && n> 1) {
		score = get_path(aln,smn,m-1,n-1,&len);
		for (i=len; i>0; i--) 
		{ int	a = aln[0][i],
			b = aln[1][i];
			sim[b][a] += smn[b][a] * PATH_WT;
		}
	}
	if (m<nb && n<na) {
		s = (float**)alloca(sizeof(float*)*(nb+2)); TEST(s)
		for (i=n; i<nb+2; i++) s[i-n] = smn[i]+m;
		score = get_path(aln,s,na-m,nb-n,&len);
		for (i=len; i>0; i--) 
		{ int	a = aln[0][i],
			b = aln[1][i];
			sim[b+n][a+m] += s[b][a] * PATH_WT;
		}
	}
	return score;
}

float get_path (aln,sim,na,nb,length) 
int	**aln;
float	**sim;
int	na, nb, *length;
{
float	**mat;
int	**ptr, i, j, k;
float	score, *colmax, rowmax, gap = 10.0; /* was 1 */
int	*maxcol, maxrow, maxi, maxj, len;
int	naa = na+2, nbb = nb+2, now;
	mat = (float**)alloca(sizeof(float*)*2); TEST(mat)
	for (i=0; i<2; i++) {
		 mat[i] = (float*)alloca(sizeof(float)*naa); TEST(mat[i])
	}
	ptr = (int**)alloca(sizeof(int*)*nbb); TEST(ptr)
	for (i=0; i<nbb; i++) {
		ptr[i] = (int*)alloca(sizeof(int)*naa); TEST(ptr[i])
	}
	colmax = (float*)alloca(sizeof(float)*(naa)); TEST(colmax)
	maxcol = (int*)alloca(sizeof(int)*(naa)); TEST(maxcol)
	for (i=0; i<naa; i++) {
		maxcol[i] = 0;
		colmax[i] = mat[0][i] = mat[1][i] = -1.0;
	}
	score = 0.0;
	now = 1;
	for (i=1; i<nbb; i++) {
		rowmax = -1.0;
		for (j=1; j<naa; j++)
		{ float dig, col, row, max;
		  int	cop, rop, top;
			rop = cop = top = 0;
			if (j>na || i>nb) mat[now][j] = 0.0;
				     else mat[now][j] = sim[i][j];
			dig = mat[!now][j-1];
			col = colmax[j-1] - gap;
			row = rowmax - gap;
			if (col > dig) {
				cop = i-maxcol[j-1]-1;
			} else {
				colmax[j-1] = dig;
				maxcol[j-1] = i-1;
			}
			if (row > dig) {
				rop = -(j-maxrow-1);
			} else {
				rowmax = dig;
				maxrow = j-1;
			}
			max = dig;
                        if (row > max) { max = row; top = rop; }
                        if (col > max) { max = col; top = cop; }
			mat[now][j] += max;
			ptr[i][j] = top;
			if (mat[now][j] > score) {
				score = mat[now][j];
				maxi = i;
				maxj = j;
			}
		}
		now = !now;
	}
	*length = 0;
	if (score > 0.1) *length = trace(sim,ptr,aln,0,maxi,maxj);
	return score;
}

trace (s,p,a,n,i,j) float **s; int **p, **a, n, i, j;
{
	if (s[i][j] < 0.0) return n;
	n++;
	a[0][n] = j;
	a[1][n] = i;
	if (i<=1 || j<=1) return n;
	if (p[i][j] > 0) i -= p[i][j];
	if (p[i][j] < 0) j += p[i][j];
	return trace(s,p,a,n,--i,--j);
}

protin (prot,seq,m,c,z,flip)
Pdbentry_ *prot;
Seq *seq; Tri ***m; Cel ***c; float z; int flip;
{
FILE	*pdb;
int	i, j, len;
Tri	**mat;
Cel	**cel;
	len = copyca(prot->Chains,seq,flip,z);
        mat = (Tri**)malloc(sizeof(Tri*)*(len+2));
        cel = (Cel**)malloc(sizeof(Cel*)*(len+2));
        for (i=0; i<=len+1; i++) {
                mat[i] = (Tri*)malloc(sizeof(Tri)*(len+2));
                cel[i] = (Cel*)malloc(sizeof(Cel)*(len+2));
        }
	add_cb(seq);
	set_vect(seq->ca,seq->cb,mat,cel,len);
	set_cbcb(seq->ca,seq->cb,mat,len);
	*m = mat;
	*c = cel;
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

celsort (ac,bc) const void *ac, *bc;
{
Cel    *a = (Cel*)ac, *b = (Cel*)bc;
        if (a->c < b->c) return 1;
        if (a->c > b->c) return -1;
        return 0;
}

set_vect (a,b,m,c,l) Vec *a, *b; Tri **m; Cel **c; int l; {
int     i, j;
Mat     frame;
float	size = 0.3, cut = 20.0;
        for (i=1; i<=l; i++) {
                setframe(a[i-1],a[i],a[i+1],&frame);
                for (j=1; j<=l; j++) { Vec s, t;
                        m[i][j].d = vdif(a[i],a[j]);
                        vinit(&(m[i][j].v));
                        if (i==j) continue;
                        vsub(a[j],a[i],&s);
                        VmulM(&frame,s,&(m[i][j].v));
			c[i][j].c = -1;
			if (fabs(m[i][j].v.x) > cut) continue;
			if (fabs(m[i][j].v.y) > cut) continue;
			if (fabs(m[i][j].v.z) > cut) continue;
			c[i][j].c = 1000000*(int)(100.0+m[i][j].v.x*size)
				   + 1000*(int)(100.0+m[i][j].v.y*size)
				     +    (int)(100.0+m[i][j].v.z*size);
			c[i][j].r = j;
		}
		c[i][0].c = c[i][j].c = -1;
		qsort(c[i],l,sizeof(Cel),celsort);
	}
}

set_cbcb (a,b,m,l) Vec *a, *b; Tri **m; int l; {
int	i, j;
	for (i=1; i<=l; i++)
	{ Vec ai, bi, ci;
		vsub(a[i+1],a[i-1],&ai);
		vnorm(&ai);
		vsub(b[i],a[i],&bi);
		vnorm(&bi);
		vprod(ai,bi,&ci);
		for (j=1; j<=l; j++)
		{ Vec aj, bj, cj;
			vsub(a[j+1],a[j-1],&aj);
			vnorm(&aj);
			vsub(b[j],a[j],&bj);
			vnorm(&bj);
			vprod(aj,bj,&cj);
			m[i][j].cos.x = vdot(ai,aj);
			m[i][j].cos.y = vdot(bi,bj);
			m[i][j].cos.z = vdot(ci,cj);
		}
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
	Vec	*ca, *cb;
	float	*acc;
	n = pdb->Aano;
	seq = (char*)malloc(sizeof(char)*(n+3));
	acc = (float*)malloc(sizeof(float)*(n+3));
	ca = (Vec*)malloc(sizeof(Vec)*(n+3));
	cb = (Vec*)malloc(sizeof(Vec)*(n+3));
	for (i=0; i<n; i++) {
		ca[i+1].x = pdb->Atoms[i].X;
		ca[i+1].y = pdb->Atoms[i].Y;
		ca[i+1].z = pdb->Atoms[i].Z;
		acc[i+1] = pdb->Atoms[i].Bfact;
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
	s->ca = ca;
	s->cb = cb;
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

moment (mom,struc,weight,natom)
float   **mom, *weight;
Vec    *struc;
int     natom;
{
        Vec   cog;
	float sum;
        int   i, j;
        vinit(&cog);
        for (i=1; i<=natom; i++)
	{ Vec	a;
		vcopy(struc[i],&a);
		vmul(&a,weight[i]);
                vsum(a,&cog);
		sum += weight[i];
        }
        vdiv(&cog,sum);
	mom[0][0] = sum;
	mom[0][1] = cog.x;
	mom[0][2] = cog.y;
	mom[0][3] = cog.z;
        for (i=1; i<4; i++) for (j=1; j<4; j++) mom[i][j] = 0.0;
        for (i=1; i<=natom; i++)
        { Vec   b;
                vsub(struc[i],cog,&b);
		vmul(&b,weight[i]);
                mom[1][1] += b.x * b.x;
                mom[2][2] += b.y * b.y;
                mom[3][3] += b.z * b.z;
                mom[1][2] = mom[2][1] += b.x * b.y;
                mom[1][3] = mom[3][1] += b.x * b.z;
                mom[2][3] = mom[3][2] += b.y * b.z;
        }
}

float super (stra, strb, seqa, seqb, sim, aln, len)
Tri	**stra, **strb;
Seq	*seqa, *seqb;
float	**sim;
int	**aln, len;
{
FILE	*out;
Sqmat_	rot;
Mat	mat;
Vec	acnt, bcnt, *axis;
int	*mm;
float	*w, *sa, *sb;
double	*ww, *ac, *bc;
double	**va, **vb;
int	na = seqa->len;
int	nb = seqb->len;
float	cut, rms, pct, suw, sud, sum = 0.0;
float	sumrms = 0.0;
int	i, j, n, id;
	sa = (float*)alloca(sizeof(float)*(na+2));
	for (i=1; i<=na; i++) sa[i] = 0.0;
	sb = (float*)alloca(sizeof(float)*(nb+2));
	for (i=1; i<=nb; i++) sb[i] = 0.0;
	mm = (int*)alloca(sizeof(int)*len);
	w = (float*)alloca(sizeof(float)*len);
	ww = (double*)alloca(sizeof(double)*len);
	va = (double**)alloca(sizeof(double*)*len);
	vb = (double**)alloca(sizeof(double*)*len);
	axis = (Vec*)alloca(sizeof(Vec)*(len+1));
	id = 0;
Pi(len) NL
	for (i=len; i>0; i--) 
	{ int	a = aln[0][i], b = aln[1][i],
		h = i-1;
	  char  ra = seqa->res[a],
		rb = seqb->res[b],
		aa1, aa2, ab1, ab2;
	  float s;
		if (ra==rb) id++;
		aa1 = aa2 = ab1 = ab2 = ' ';
		if (seqa->acc[a]>0.0) aa1 = '*';
		if (seqb->acc[b]>0.0) ab1 = '*';
		if (seqa->acc[a]>0.5) aa2 = '*';
		if (seqb->acc[b]>0.5) ab2 = '*';
		s = 1.0;//rescore(stra,strb,a,b,aln,len); 
		sum += s;
		printf("%c%c%c %4d %5.1f%4d %c%c%c\n",
			aa2,aa1,ra,a,s,b,rb,ab1,ab2);
		va[h] = (double*)alloca(sizeof(double)*3);
	 	va[h][0] = seqa->ca[a].x;
	 	va[h][1] = seqa->ca[a].y;
	 	va[h][2] = seqa->ca[a].z;
		vb[h] = (double*)alloca(sizeof(double)*3);
	 	vb[h][0] = seqb->ca[b].x;
	 	vb[h][1] = seqb->ca[b].y;
	 	vb[h][2] = seqb->ca[b].z;
		sa[a] = sb[b] = sqrt(s);
		w[i] = s;
	}
	pct = 100.0*(float)id/(float)len;
	for (i=0; i<10; i++) norm(3.0,w+1,len);
	suw = 0.0;
	for (i=0; i<len; i++) {
		w[i] = ww[i] = (double)w[i+1];
		suw += w[i];
	}
	rot = alloc_sqmat(3);
	ac = (double*)alloca(sizeof(double)*3);
	bc = (double*)alloca(sizeof(double)*3);
	rms = supermac(ww,va,vb,len,ac,bc,rot);
	acnt.x = ac[0]; acnt.y = ac[1]; acnt.z = ac[2]; 
	bcnt.x = bc[0]; bcnt.y = bc[1]; bcnt.z = bc[2]; 
	mat.A.x = rot[0][0]; mat.A.y = rot[0][1]; mat.A.z = rot[0][2];
	mat.B.x = rot[1][0]; mat.B.y = rot[1][1]; mat.B.z = rot[1][2];
	mat.C.x = rot[2][0]; mat.C.y = rot[2][1]; mat.C.z = rot[2][2];
	for (i=0; i<=na+1; i++) vsub(seqa->ca[i],acnt,seqa->ca+i);
	for (i=0; i<=nb+1; i++) {
		vsub(seqb->ca[i],bcnt,seqb->ca+i);
		MmulV(&mat,seqb->ca[i],seqb->ca+i);
	}
	sud = 0.0;
	for (i=len; i>0; i--) 
	{ int	a = aln[0][i], b = aln[1][i];
	  float d = vddif(seqa->ca[a],seqb->ca[b]);
		sud += exp(-d*0.1);
	}	
	sud = 100*sud/(float)len;
	printf("%s Weighted RMSd = %7.3f over %d atoms %f %f %5.1f %f\n",
		pdbcode, rms, len, suw, sum, pct, sud);
/*
	if (self)
	{ float wmax, **mom, **vec, val[4], x, y, z;
	  Vec	cog;
	  FILE  *axs;
		axs = fopen("axis.pdb","w");
		n = 0;
		for (i=len; i>0; i--) 
		{ int	a = aln[0][i], b = aln[1][i];
			n++;
			vave(seqa->ca[a],seqa->ca[b],axis+i);
             		fprintf(axs,"ATOM%7d  CA  GLY C%4d     %7.3f %7.3f %7.3f  0.00 %5.2f\n",
                        	n, n, axis[i].x, axis[i].y, axis[i].z, w[i]);
		}
                fprintf(axs,"TER\n");
                mom = (float**)alloca(sizeof(float*)*4);
                vec = (float**)alloca(sizeof(float*)*4);
        	for (i=0; i<4; i++) {
                	mom[i] = (float*)alloca(sizeof(float)*4);
                	vec[i] = (float*)alloca(sizeof(float)*4);
		}
		moment(mom,axis,w,len);
		eigen(mom,3,val,vec);
		x = mom[0][1]-vec[1][1]*10.0;
		y = mom[0][2]-vec[2][1]*10.0;
		z = mom[0][3]-vec[3][1]*10.0;
		n = 1;
                fprintf(axs,"ATOM%7d  CA  GLY D%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", n, n, x, y, z);
		cog.x = x = mom[0][1];
		cog.y = y = mom[0][2];
		cog.z = z = mom[0][3];
		n = 2;
                fprintf(axs,"ATOM%7d  CA  GLY D%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", n, n, x, y, z);
		x = mom[0][1]+vec[1][1]*10.0;
		y = mom[0][2]+vec[2][1]*10.0;
		z = mom[0][3]+vec[3][1]*10.0;
		n = 3;
                fprintf(axs,"ATOM%7d  CA  GLY D%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", n, n, x, y, z);
                fprintf(axs,"TER\n");
		fclose(axs);
	}
*/
	superout(seqa, seqb, sa, sb, aln, len);
	cut = 0.5;
	for (j=15; j; j--) {
		n = 0;
		for (i=0; i<len; i++) {
			if (w[i] > cut) {
				n++; ww[i] = 1.0;
			} else  ww[i] = 0.0;
		}
		rms = supermac(ww,va,vb,len,ac,bc,rot);
		printf("%s Un-weighted RMSd = %7.3f over best %d atoms ( %f  %5.1f )\n",
			pdbcode, rms, n, sum, pct);
		cut *= 0.5;
	}
	out = fopen("plot.rms","w");
	sort(0,w,0,mm,len,1);
	sumrms = 0.0;
	for (i=5; i<len; i++) {
		for (j=0; j<len; j++) ww[j] = 0.0;
		for (j=0; j<=i; j++) ww[mm[j]] = 1.0;
		rms = supermac(ww,va,vb,len,ac,bc,rot);
		sumrms += rms;
		fprintf(out,"%d %7.3f %7.3f\n", j, rms, sumrms);
	}
	n = 0;
	for (i=0; i<len; i++) {
		n++;
		ww[i] = 1.0;
	}
	rms = supermac(ww,va,vb,len,ac,bc,rot);
	printf("%s Un-weighted RMSd = %7.3f over all %d matched atoms %f %f  %5.1f %f\n",
		pdbcode, rms, n, suw, sum, pct, sud);
	return rms;
}

superout (seqa, seqb, sa, sb, aln, len)
Seq	*seqa, *seqb;
float	*sa, *sb;
int	**aln, len;
{
FILE	*out;
int     i, nb, n = 1;
int     lena, lenb;
char    aa3[80], aaa[4];
        strcpy(aa3,
	"ALAASXCYSASPGLUPHEGLYHISILEACELYSLEUMETASNPCAPROGLNARGSERTHRUNKVALTRPXXXTYRGLX");
	out = fopen("super.pdb","w");
	lena = seqa->len;
        for (i=1; i<=lena; i++) {
		strncpy(aaa,aa3+3*(seqa->res[i]-'A'),3); aaa[3] = 0;
                fprintf(out,"ATOM%7d  CA  %s A%4d     %7.3f %7.3f %7.3f %5.2f %5.2f\n",
                        n, aaa, i, seqa->ca[i].x, seqa->ca[i].y, seqa->ca[i].z,
			0.5*seqa->acc[i]+0.5, sa[i]);
		n++;
        }
        fprintf(out,"TER\n");
	nb = n;
	lenb = seqb->len;
        for (i=1; i<=lenb; i++) {
		strncpy(aaa,aa3+3*(seqb->res[i]-'A'),3); aaa[3] = 0;
                fprintf(out,"ATOM%7d  CA  %s B%4d     %7.3f %7.3f %7.3f %5.2f %5.2f\n",
                        n, aaa, i, seqb->ca[i].x, seqb->ca[i].y, seqb->ca[i].z,
			0.5*seqb->acc[i]+0.5, sb[i]);
		n++;
        }
        fprintf(out,"TER\n");
	for (i=len; i>0; i--) 
	{ int	a = aln[0][i],
		b = aln[1][i] + nb - 1;
                fprintf(out,"CONECT  %3d    0    0    0    0  %3d", a,b);
                fprintf(out,"                                           :\n");
	}
        fprintf(out,"END\n");
	fclose(out);
}

/*
matplot("bias.ps",bias,lenb,lena,0);
exit(1);
*/
matplot (file,mat,lena,lenb,logs) 
char	*file;
float	**mat;
int	lena,lenb, logs;
{
float	hi, lo, span;
int	i, j, size=3;
FILE	*out;
	hi = -99999.9; lo =  99999.9;
	for (i=0; i<lena; i++) {
		for (j=0; j<lenb; j++) { float x = mat[i][j];
			if (x>hi) hi = x;
			if (x<lo) lo = x;
		}
	}
	span = hi-lo;
	out = fopen(file,"w");
	fprintf(out,"/Times-Roman findfont\n12 scalefont\nsetfont\n");
	fprintf(out,"newpath\n");
	for (i=0; i<lena; i++) {
		for (j=0; j<lenb; j++) { float r,g,b, x;
			x = (mat[i][j]-lo)/span;
			if (logs) x = 3.0*log(x+1.0)/log(2)-1.0;
				else x = 3.0*x-1.0;
			if (x>0) {
				if (x>1.0) { r = x-1.0; g = 2.0-x; b = 0.0; }
			      		else { b = 1.0-x; g = x; r = 0.0; }
			} else { r = 0.0; g = 0.0; b = 1.0+x; }
			fprintf(out,"%d %d moveto\n", i*size,j*size);
			fprintf(out,"0 10 rlineto\n10 0 rlineto\n0 -10 rlineto\nclosepath\n");
			fprintf(out,"%6.3f %6.3f %6.3f setrgbcolor\nfill\n", r,g,b);
		}
	}
	fprintf(out,"showpage\n");
	fclose(out);
}
