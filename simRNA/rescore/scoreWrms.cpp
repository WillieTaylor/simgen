/*
c++ scoreWrms.cpp -o scoreWrms sims/util.o sims/geom.o -lm

./scoreWrms model.cas psicov.top (needs true.pdb)

*/
#include "sims/util.hpp"
#include "sims/geom.hpp"

#include "util/aa/incl/matrix.h"
#include "util/aa/incl/bestrot.h"
#include "util/aa/incl/siva.h"
#include "util/aa/incl/ql.h"

typedef struct {
	char	*res;	// residue type
	float	*acc;	// surface accesssibility
	int	*sec;	// secondary structure type (0=c, 1=A, 2=B, 3=3ten)
	int	*rid;	// residue number
	Vec	*cas;	// CA position (runs 1...N with dummy 0 and N+1)
	Vec	*cbs;	// dummy CB/centroid (extended <bond>A from CA)
	int	len;	// chain length
	int	gaps;	// number of chain breaks
	float	pcta, pctb;	// percent alpha and beta structure
} Prot;

Vec	rms121 ( Prot*, Prot* );
void	smooth ( Vec*, int, int );
double	supermac ( double*, double**, double**, int, double*, double*, double** );

char	aa3[2222][4];

void addres ( Vec* cas, int i, int j, int k, int add ) {
// Add a new residue position <add> off the end of the chain away from mid <j> and <k>
Vec	m = cas[j] & cas[k];
Vec	v = m - cas[i];
	cas[add] = m+v;
}

void extend ( Prot* chain ) {
// Add dummy CA atoms to the end of the chain
int	n = chain->len, m = n+1; 
        addres(chain->cas,  3,  2,1, 0);	// add res 0 to Nterm
        addres(chain->cas,n-2,n-1,n, m);	// add res N+1 to Cterm
	chain->res[0] = 'n';
	chain->res[m] = 'c';
}

void add_cb ( Prot *chain, float bond ) {
// Add dummy CB atoms to the chain (but not to 0 or N+1)
        DO1(i,chain->len) { Vec n,c,b; float d;
		n = (chain->cas[i] - chain->cas[i-1]).norm();
		c = (chain->cas[i] - chain->cas[i+1]).norm();
		b = (n+c).norm()*bond;
		chain->cbs[i] = chain->cas[i]+b;
        }
}

void getpdb ( Prot* chain, FILE *pdb ) { 
// Read a PDB format file from stream <pdb> into structure Prot 1..N
int	io, n = 1;
char	line[222], junk[30];
Vec	xyz[2222];
float	acc[2222];
        while(1) { float x,y,z, s,t; int a,b; char c;
		io = read_line(pdb,line);
		if ( io < 1 ) break;
		if (!strncmp(line,"END",3)) break;
		if (!strncmp(line,"TER",3)) continue; // LINKs may follow
		if (strncmp(line,"ATOM   ",7) && strncmp(line,"HETATM ",7)) continue; // MSE is a HETATM
		if (strncmp(line+13,"CA ",3)) continue;
		sscanf(line,"%30c %f%f%f%f%f", junk, &x, &y, &z, &s, &t);
		xyz[n].x = x; xyz[n].y = y; xyz[n].z = z; acc[n] = t;
		strncpy(aa3[n],junk+17,3);
		n++;
	}
	chain->len = n-1;
	n++; // 2 extra positions for dummy chain extensions
	chain->cas = new Vec[n];
	chain->cbs = new Vec[n];
	chain->sec = new int[n];
	chain->rid = new int[n];
	chain->res = new char[n];
	chain->acc = new float[n];
	DO1(i,chain->len) {
		chain->cas[i].x = xyz[i].x;
		chain->cas[i].y = xyz[i].y;
		chain->cas[i].z = xyz[i].z;
		chain->acc[i] = acc[i];
	}
	extend(chain);		// add 0 and N+1 positions
	add_cb(chain,2.0);	// add pseudo CB/centroid
}

float** getmat ( Prot *chain )
{
int	len = chain->len, lenn = len+2; // allow for dummy atoms 0 and N+1
float	**mat = new float*[lenn];
	DO(i,lenn) mat[i] = new float[lenn];
	DO(i,lenn) DO(j,lenn) {
		mat[i][j] = 0.0;
		if (i<j) mat[i][j] = chain->cas[i] | chain->cas[j];
		if (i>j) mat[i][j] = chain->cbs[i] | chain->cbs[j];
	}
	return mat;
}

Vec rogs ( Prot *A, float *w ) {
// get the weighted RoG around a rough bundle axis
int	m, n, in, len = A->len;
Vec	*c = new Vec[len];
Pairs	*p = new Pairs[999];	
int	q[999];
Seg	axis;
Vec	cog, move;
Vec	sum, rog;
	m = n = 0;
	cog.zero();
	DO1(i,len) { // caps have weight = 1
		if ( w[i]>0.5) { cog += A->cas[i]; m++; }
		if ( w[i]>0.5 && w[i]<1.5) c[n++] = A->cas[i];
	}
	cog /= (float)m;
	m = 0;
	DO(i,n-8) { // for pairs of cap residues...
		for (int j=i-1; j<i+8; j++) // between adjacent caps
		{ float	d = c[i] | c[j];
			if (d < 30.0) continue;
			p[m].a = i; p[m].b = j; p[m].s = d; m++;
			if (m==999) break;
		}
	}
	sort(p,q,m); // reverse sort with rank in <q[]>
	axis.A = c[p[q[0]].a]; axis.B = c[p[q[0]].b]; // first axis guess is widest cap pair
	in = 1;
	DO(j,m)
	{ int	i = q[j], pa = p[i].a, pb = p[i].b;
	  Vec	a = c[pa], b = c[pb];
	  Seg	now = (axis.A/(float)j, axis.B/(float)j);
	  float dA = a|now.A, dB = a|now.B;
		if (dA < dB) {
			axis.A += a; axis.B += b; 
		} else {
			axis.A += b; axis.B += a;
		} 
	}
	axis.A /= (float)m;
	axis.B /= (float)m;
	move = cog-(axis.A & axis.B);
	axis.A += move; axis.B += move;
	// the bundle axis is now centred on the CoG of the TM segments
	rog.zero(); sum.zero();
	DO1(i,len)
	{ float d = axis.vec_to_line(A->cas[i]),
		dd = d*d, wi = w[i];
		rog.x += dd;    sum.x += 1.0;
		rog.y += dd*wi; sum.y += wi;
		wi += 0.5;
		rog.z += dd*wi; sum.z += wi;
	}
	rog.x /= sum.x; rog.x = sqrt(rog.x);
	rog.y /= sum.y; rog.y = sqrt(rog.y);
	rog.z /= sum.z; rog.z = sqrt(rog.z);
	return rog;
}

// 1 = pdb file, 2 = constraints file

int main (int argc, char** argv) {
// XYZ from simprot realigned with tube axis along Z (1st pdb line = prot.dat top GROUP line)
Vec	rms, rog;
float	sum0, sum1, sum2;
int	n, len, in=500;
char	line[111];
Pairs	cons[1111];
float	**dij;
Prot	*chain = new Prot;
Prot	*known = new Prot;
FILE	*con = fopen(argv[1],"r");
FILE	*pdb = fopen(argv[2],"r");
float	want, give, bump = 4.0;
int	nbump = 0;
	getpdb(chain,pdb);
	fclose(pdb);
	pdb = fopen("true.pdb","r");
	getpdb(known,pdb);
	fclose(pdb);
	dij = getmat(chain); // i<j=CA i>j=CB
	len = chain->len;
	DO(i,len) DO(j,len) { float d;
		if (i<j) continue;
		if (i-j < 5) continue;
		d = chain->cas[i] | chain->cas[j];
		if (d < bump) nbump++;
	}
	Pi(nbump) NL
	n = 0;
	DO(i,in) { int a,b,c,d; float s;
		if (read_line(con,line) <= 0) break;
		sscanf(line,"%d %d %d %d", &a, &b, &c, &d);
		if (a>b) { int e=a; a=b; b=e; } // keep a<b
		if (d>9000) d -= 9000; // added to stop early breaks in sim
		s = 0.001*(float)d;
		cons[n].a = a;
		cons[n].b = b;
		cons[n].s = s;
		n++;
	}
	sscanf(argv[3],"%f", &want);
	sscanf(argv[4],"%f", &give);
	printf("Ideal distance = %f, +/-%f\n", want, give);
	give = 1.0/(give*give);
	sum0 = sum1 = sum2 = 0.0;
	DO(i,n) { int a = cons[i].a, b = cons[i].b; float s = cons[i].s, d, e, wa,wb;
		d = dij[a][b]-want; e = exp(-d*d*give); // CB=[b][a], CA=[a][b]
		//Pi(i) Pi(a) Pi(b) Pr(s) Pr(dij[a][b]) Pr(want) Pr(d) Pr(e) NL
		// Pi(i) Pi(a) Pi(b) Pr(s) Pr(dij[a][b]) Pr(dij[b][a]) Pr(e) NL
		wa = chain->acc[a]; wb = chain->acc[b];
Pi(a) Pi(b) Pr(d) Pr(e) Pr(s) Pr(wa) Pr(wb) 
		sum0 += s*e;
		sum1 += s*e*wa*wb;
		sum2 += s*e*(wa+0.5)*(wb+0.5);
		Pi(i+1) Pr(sum0) Pr(sum1) Pr(sum2) NL
	}
	rms = rms121(chain,known);
	printf("RMS   Unsmoothed %6.3f  Smoothed 3x %6.3f  Smoothed 6x %6.3f\n", rms.x, rms.y, rms.z);
	pdb = fopen("smoothA.pdb","w");
	DO1(i,chain->len) {
		fprintf(pdb,"ATOM%7d  CA  GLY A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,
			chain->cas[i].x, chain->cas[i].y, chain->cas[i].z);
	}
	fclose(pdb);
	pdb = fopen("smoothB.pdb","w");
	DO1(i,known->len) {
		fprintf(pdb,"ATOM%7d  CA  GLY A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,
			known->cas[i].x, known->cas[i].y, known->cas[i].z);
	}
	fclose(pdb);
}

Vec rms121 ( Prot *a, Prot *b )
{
Vec	rms;
double  **rot;
double	*ww, *ac, *bc;
double	**va, **vb;
int	i, j, len = a->len;
int	cycles = 3;
	rot = (double**)alloca(sizeof(double*)*3);
	DO(i,3) rot[i] = (double*)alloca(sizeof(double)*3);
	ac = (double*)alloca(sizeof(double)*3);
	bc = (double*)alloca(sizeof(double)*3);
	ww = (double*)alloca(sizeof(double)*len);
	va = (double**)alloca(sizeof(double*)*len);
	vb = (double**)alloca(sizeof(double*)*len);
	DO1(i,len) { int h = i-1;
		va[h] = (double*)alloca(sizeof(double)*3);
		vb[h] = (double*)alloca(sizeof(double)*3);
	}
	DO1(i,len) { int h = i-1;
	 	va[h][0] = a->cas[i].x; va[h][1] = a->cas[i].y; va[h][2] = a->cas[i].z;
	 	vb[h][0] = b->cas[i].x; vb[h][1] = b->cas[i].y; vb[h][2] = b->cas[i].z;
	}
	// Unsmoothed
	DO(i,len) ww[i] = b->acc[i];
	rms.x = supermac(ww,va,vb,len,ac,bc,rot);
	if (rms.x>99.999) rms.x = 99.999;
	// Smoothed 3x
	smooth(a->cas,len,cycles);
	smooth(b->cas,len,cycles);
	DO1(i,len) { int h = i-1;
	 	va[h][0] = a->cas[i].x; va[h][1] = a->cas[i].y; va[h][2] = a->cas[i].z;
	 	vb[h][0] = b->cas[i].x; vb[h][1] = b->cas[i].y; vb[h][2] = b->cas[i].z;
	}
	rms.y = supermac(ww,va,vb,len,ac,bc,rot);
	if (rms.y>99.999) rms.y = 99.999;
	// Smoothed 6x
	smooth(a->cas,len,cycles);
	smooth(b->cas,len,cycles);
	DO1(i,len) { int h = i-1;
	 	va[h][0] = a->cas[i].x; va[h][1] = a->cas[i].y; va[h][2] = a->cas[i].z;
	 	vb[h][0] = b->cas[i].x; vb[h][1] = b->cas[i].y; vb[h][2] = b->cas[i].z;
	}
	rms.z = supermac(ww,va,vb,len,ac,bc,rot);
	if (rms.z>99.999) rms.z = 99.999;
	return rms;
}

void smooth ( Vec *a, int len, int cycles)
{
	DO(i,cycles)
	{ Vec	p1 = a[0], p2 = a[1], p3 = a[2];
		DO1(j,len) { Vec q;
			q = (p1+p2+p3)/3;
			p1 = a[j]; p2 = a[j+1]; p3 = a[j+2];
			a[j] = q;
		}
	}
}


/* ==== FUNCTIONS bestrot.c ==== */

/* An implementation of the 3D point set alignment algorithm by
 * A. D. McLachlan. Reference:
 * McLachlan, A. D. (1979): J. Mol. Biol. 128: 49-79.
 * Replaces the buggy Kabsch rotation algorithm.
 */

/* ANSI C, IRIX 5.2, 5. Aug. 1994. Andris Aszodi */

/* ---- HEADER ---- */

//#include "incl/bestrot.h"

/* ---- INCLUDE FILES ---- */

//#include "incl/siva.h"	    /* singular value decomposition */

/* ---- DEFINITIONS ---- */

#define DIM 3

double supermac (double *W, double **X, double **Y, int n,
		 double *CtrX, double *CtrY, double **rot)
//		 double *CtrX, double *CtrY, Sqmat_ rot)
{
double	rms;
	center_vectors(X,CtrX,W,n);
	center_vectors(Y,CtrY,W,n);
	rms = best_rot(X,Y,W,n,rot);
	if (rms < 0.0) {
		X[0][0] += drand48();
		X[0][1] += drand48();
		X[0][2] += drand48();
		Y[0][0] += drand48();
		Y[0][1] += drand48();
		Y[0][2] += drand48();
		rms = best_rot(X,Y,W,n,rot);
		printf("*NB* tried to patch-up bad super (rms=%f)\n", rms);
	}
	return rms;
}

/* ==== FUNCTIONS ==== */

/* center_vectors: calculates the centroid of the 
 * set of 3-dimensional vectors X (Vno x 3) and
 * subtracts it from each of them,  thus centring the
 * set on the centroid. If Ctr==NULL, then a 3-long
 * array is allocated to store the centroid coordinates;
 * if Ctr!=NULL, then it is assumed to be large enough to
 * hold the coordinates.
 * Return value: Ctr, or NULL if Vno==0.
 */
double *center_vectors(double **X, double *Ctr, const double *W, unsigned int Vno)
{
    register unsigned int i, j;
    double	Wsum;
    
    if (!Vno) return(NULL);
    
    /* allocate centroid vector if absent */
    if (Ctr==NULL)
	Ctr=(double *) calloc(DIM, sizeof(double));
    
    /* get the (weighted) centroid */
    Wsum = 0.0;
    for (i=0; i<Vno; i++) Wsum+=W[i];
    for (j=0; j<DIM; j++)
    {
	Ctr[j]=0.0;
	for (i=0; i<Vno; i++) Ctr[j]+=X[i][j]*W[i];
	Ctr[j]/=Wsum;
    }
    
    /* subtract Ctr from all vectors in X */
    for (i=0; i<Vno; i++)
	for (j=0; j<DIM; j++)
	    X[i][j]-=Ctr[j];
    
    return(Ctr);
}
/* END of center_vectors */

/* best_rot: finds the best rotation matrix that brings a set of
 * vectors X into another set Y. X, Y have Vno vectors (in rows), 
 * and both live in 3 dimensions (Vno x 3). W is a Vno-long
 * weight vector that can emphasise vector pairs. Transform is
 * a 3x3 square matrix (allocated before call) that on
 * return contains the X->Y transformation. It is assumed that
 * X and Y were centered before the call.
 * NOTE: this routine cannot handle the degenerate cases when
 * the point sets are Dim<3-dimensional. (Might be implemented
 * later.) When this happens, a warning is printed to stderr
 * and -1.0 (a meaningless RMS value) is returned.
 * Return value: a weighted least-squares error function.
 */
double best_rot(double **X, double **Y, const double *W, 
		unsigned int Vno, Sqmat_ Transform)
{
    register unsigned int i, j, k, n, m;
    int Psign, Rank;
    double **H=NULL, *D=NULL, **K=NULL;
    Sqmat_ U=NULL;
    register double Err=0.0, Temp1, Temp2, Wsum=0.0;
    double Detu;
    
    /* set up the matrix to be SVD-d */
    U=alloc_sqmat(DIM);
    for (i=0; i<DIM; i++)
	for (j=0; j<DIM; j++)
	{
	    Temp1=0.0;
	    for (k=0; k<Vno; k++) Temp1+=W[k]*X[k][i]*Y[k][j];
	    U[i][j]=Temp1;
	}
    
    /* set up and perform SVD */
    siva_setup(DIM, DIM, &H, &D, &K);
    siva_decomp((const double**)U, DIM, DIM, H, D, K);
    
    /* check rank: do nothing if rank was lost */
    if ((Rank=rank_cond(D, DIM, SIVA_EPSILON, NULL))<3)
    {
	fprintf(stderr, "? best_rot(): Rank %d<%d\n", Rank, DIM);
	free_siva(DIM, DIM, H, D, K); free_matrix(U, DIM);
	return(-1.0);
    }
    
    /* get the determinant of U and store its sign in Psign */
    Psign=lu_decomp(U, DIM, NULL);
    Detu=lu_det(U, Psign, DIM);
    Psign=(Detu>0)? 1: -1;
    
    /* generate the transform matrix: here we explicitly
     * use DIM==3 because McLachlan does not say what to do
     * if Psign==-1 and DIM>3 and I don't know :-)
     */
    for (i=0; i<DIM; i++)
	for (j=0; j<DIM; j++)
	    Transform[i][j]=
		K[i][0]*H[j][0]+K[i][1]*H[j][1]+Psign*K[i][2]*H[j][2];
    
    /* evaluate the error function */
    for (n=0; n<Vno; n++)
    {
	Temp2=0.0;
	for (i=0; i<DIM; i++)
	{
	    Temp1=0.0;
	    for (j=0; j<DIM; j++)
		Temp1+=Transform[i][j]*X[n][j];
	    Temp1-=Y[n][i];
	    Temp2+=Temp1*Temp1;
	}
	Err+=W[n]*Temp2;
	Wsum+=W[n];
    }
    Err/=Wsum;
    
    /* cleanup */
    free_siva(DIM, DIM, H, D, K); free_matrix(U, DIM);
    return(sqrt(Err));
}
/* END of best_rot */

#undef DIM

/* ==== END OF FUNCTIONS bestrot.c ==== */

/* ==== FUNCTIONS siva.c ==== */

/* Singular value decomposition: homebrew version */

/* ANSI C, IRIX 4.0.5, 21. July 1994. Andris Aszodi */

/* ----	HEADER AND INCLUDE FILES ---- */

//#include "incl/matrix.h"
//#include "incl/ql.h"
//#include "incl/siva.h"

/* ==== FUNCTIONS ==== */

/* ---- SINGULAR VALUE DECOMPOSITION AND BACK-SUBSTITUTION ---- */

/* siva_setup: allocates the matrices etc. for storing the
 * decomposed form. The matrix to be decomposed is (Row x Col), 
 * where Row>=Col. If Row<Col, then the rows in U will be padded
 * to get a Col x Col matrix. The actual Row is returned.
 * U is (Row x Col), W is Col long, V is (Col x Col).
 */
int siva_setup(int Row, int Col, double ***U, double **W, double ***V)
{
    double **u, *w, **v;    /* local matrices and vectors */
    int i;
    
    if (Row<=0 || Col<=0)	/* para */
    {
	fprintf(stderr, "? siva_setup(): Bad dimensions (%d x %d)\n", 
		Row, Col);
	return(Row);
    }
    
    /* If Row<Col, then u is padded with Col-Row zeroed rows 
     * to make it a Col x Col matrix 
     * and a warning is printed. 
     */
    if (Row<Col)
    {
	fprintf(stderr, "? siva_setup(): %d x %d matrix, rows padded\n", 
	    Row, Col);
	Row=Col; 
    }
    u=(double **) calloc(Row, sizeof(double *));
    for (i=0; i<Row; i++)
	u[i]=(double *) calloc(Col, sizeof(double));
    
    w=(double *) calloc(Col, sizeof(double));
    v=(double **) calloc(Col, sizeof(double *));
    for (i=0; i<Col; i++)
	v[i]=(double *) calloc(Col, sizeof(double));
    
    /* return */
    *U=u; *W=w; *V=v;
    return(Row);    /* actual no. of rows, may be padded */
}
/* END of siva_setup */

/* siva_decomp: singular value decomposition routine 
 * based on a Hungarian linear algebra book by Pa'l Ro'zsa.
 * It is assumed that
 * siva_setup() has been called prior to siva_decomp(). The "padding"
 * is taken care of automatically so Row<Col entries are accepted.
 * (DO NOT use the Row value returned by siva_setup()! )
 * SVD is carried out on A (which will be preserved).
 * The result of the SVD A=UWV' 
 * is returned in the arrays allocated by siva_setup():
 * arrays: U is a Row x Col (or Col x Col) matrix, W contains the
 * singular values in a vector and V is a Col x Col matrix.
 * Return value: 1 if the built-in iteration limit in eigen_ql()
 * is exceeded, 0 if OK.
 */
int siva_decomp(const double **A, int Row, int Col, 
	double **U, double *W, double **V)
{
    Trimat_ Ata=NULL;
    register double Temp;
    register int i, j, k;
    int Err=0;
    
    /* construct the A'A matrix lower triangle in Ata */
    Ata=alloc_trimat(Col);
    for (i=0; i<Col; i++)
	for (j=0; j<=i; j++)
	{
	    Temp=0.0;
	    for (k=0; k<Row; k++) Temp+=A[k][i]*A[k][j];
	    Ata[i][j]=Temp;
	}
     
    /* get the eigenvalues and eigenvectors of Ata:
    * W holds the eigenvalues sorted in decreasing
    * order, V holds the corresponding eigenvectors
    * as ROWS (also sorted)
    */
    Err=eigen_ql(Ata, Col, W, V);
    free_matrix(Ata, Col);
    if (Err) return(1);
    
    /* transpose V in place and sqrt the singular value vector W */
    for (i=0; i<Col; i++)
    {
	Temp=W[i];
	W[i]=(Temp<SIVA_MINVAL)? 0.0: sqrt(Temp);
	for (j=0; j<i; j++)
	{ Temp=V[i][j]; V[i][j]=V[j][i]; V[j][i]=Temp; }
    }
    
    /* get the matrix U: Ro'zsa says that A*v(j)=W[j]*u(j),
    * where u(j) and v(j) are the j-th columns of U and V, 
    * respectively. If W[j] is too small, then the corresponding
    * u(j) will be set to 0 (cf. sqrt(W) above)
    */
    for (j=0; j<Col; j++)
    {
	 if (W[j]==0.0)	/* set U[][j] to zero */
	 {
	     for (i=0; i<Row; i++) U[i][j]=0.0;
	     continue;
	 }
	 
	 for (i=0; i<Row; i++)
	 {
	     Temp=0.0;
	     for (k=0; k<Col; k++) Temp+=A[i][k]*V[k][j];
	     U[i][j]=Temp/W[j];
	 }
	 
	 /* If Row<Col, the rest of the U rows should be zero.
	  * siva_setup() does this but as a precaution we
	  * zero these rows here anyway
	  */
	 if (Row<Col)
	     for (i=Row; i<Col; i++)
		memset(U[i], 0, Col*sizeof(double));

    }
    return(0);
}
/* END of siva_decomp */

/* rank_cond: checks the N singular values W[] of a matrix 
 * after SVD. If Cond!=NULL, then the condition number 
 * (ratio of the largest and smallest singular value) is
 * calculated. The singular values which are smaller than
 * Eps times the largest are set to 0.0.
 * Return value: the rank of the matrix.
 */
int rank_cond(double W[], int N, double Eps, double *Cond)
{
    double Wmax=-HUGE_VAL, Wmin=HUGE_VAL;
    int i, Rank;
    
    /* get the largest and smallest singular value */
    for (i=0; i<N; i++)
    {
	if (W[i]>Wmax) Wmax=W[i];
	if (W[i]<Wmin) Wmin=W[i];
    }
    
    /* calc the condition number: set to HUGE_VAL if Wmin==0.0 */
    if (Cond!=NULL)
	*Cond=(Wmin==0.0)? HUGE_VAL: Wmax/Wmin;
    
    /* set all singular values which are smaller than Eps*Wmax
     * to zero: this is the conditioning. Calc the rank
     */
    Wmax*=fabs(Eps); Rank=N;
    for (i=0; i<N; i++)
	if (W[i]<Wmax) 
	{
	    W[i]=0.0; Rank--;
	}
    return(Rank);
}
/* END of rank_cond */

/* siva_solve: back-substitution routine for solving linear equations
 * AX=B. A should be SV-decomposed into U, W and V' by siva_comp()
 * and the weight vector should be "conditioned" (small entries
 * zeroed) by rank_cond() prior to the call to this routine.
 * U, W and V' are the decomposed and conditioned bits, Row and Col
 * are the row and column numbers of A (Row may have been "padded"
 * by siva_comp()!), B[] is the Row long right-hand-side vector.
 * (Padding may be necessary here, too!)
 * The result is returned in X[] (Col long).
 */
void siva_solve(const double **U, const double W[], const double **V, 
	int Row, int Col, const double B[], double X[])
{
    register double *Tmp=NULL;
    register double Sum;
    register int i, j, k;
    
    /* get WU'*B first */
    Tmp=(double *) calloc(Col, sizeof(double));	/* calloc() zeroing is essential here */
    for (j=0; j<Col; j++)
	if (W[j]!=0.0)	/* skip zeroed */
	{
	    Sum=0.0;
	    for (i=0; i<Row; i++) Sum+=U[i][j]*B[i];
	    Tmp[j]=Sum/W[j];
	}
    
    /* multiply Tmp by V to get solution vector */
    for (i=0; i<Col; i++)
    {
	Sum=0.0;
	for (j=0; j<Col; j++) Sum+=V[i][j]*Tmp[j];
	X[i]=Sum;
    }
    free(Tmp);
}
/* END of siva_solve */

/* free_siva: cleans up the space allocated to the 3 SVD arrays. */
void free_siva(int Row, int Col, double **U, double *W, double **V)
{
    int i;
    
    /* free up U: if Row<Col, then it is Col x Col */
    if (Row<Col) Row=Col;
    for (i=0; i<Row; i++) free(U[i]);
    free(U);
    for (i=0; i<Col; i++) free(V[i]);
    free(V);
    free(W);
}
/* END of free_siva */

/* ==== END OF FUNCTIONS siva.c ==== */

/* ==== FUNCTIONS matrix.c ==== */

/* Routine collection for square and lower triangle matrices: the latter
 are stored economically. */

/* ANSI C, IRIX 5.2, 5. Aug 1994. Andris Aszodi */

/* ---- HEADER ---- */

//#include "incl/matrix.h"

/* ---- DEFINITIONS ---- */

#ifdef FLT_MIN
#define LU_EPSILON (10.0*FLT_MIN)
#else
#define LU_EPSILON (1.0e-30)
#endif

/* ==== FUNCTIONS ==== */

/* ---- LOWER TRIANGLE MATRICES ---- */

/* alloc_trimat: allocates space for a triangular matrix with Size rows.
 The triangle contains the main diagonal as well. 
 Returns the pointer to the matrix or NULL if alloc failed. */

extern Trimat_ alloc_trimat(int Size)
{
    Trimat_ Mat;
    int i;
  
    Mat=(double **) calloc(Size,sizeof(double *));
    if (Mat!=NULL)
	for (i=0; i<Size; i++)
	{
	    Mat[i]=(double *)calloc(i+1,sizeof(double));
	    if (Mat[i]==NULL)
	    {	/* clean up already allocated rows */
		free_matrix(Mat,i);
		return(NULL);
	    }
	}
    return(Mat);
}
/* END of alloc_trimat */

/* free_matrix: frees up space occupied by Mat. The same routine for
 triangular and square matrices. */

extern void free_matrix(double **Mat, int Size)
{
    int i;
  
    for (i=Size-1; i>=0; i--)	/* dealloc in reverse order */
	free(Mat[i]);
    free(Mat);
}
/* END of free_matrix */

/* list_trimat: lists Mat to stdout with entries occupying Width chars,
 Prec digits precision. If a row takes up more than Linewidth chars,
 then the matrix is cut up nicely. */

extern void list_trimat(Trimat_ Mat, int Size, int Linewidth,
	int Width, int Prec)
{
    char Entrys[10],Cols[8],Rows[8];	/* format strings and len determ */
    int i,j,k,Sizew,Chunk,Jbeg,Items,Ulinelen;

    sprintf(Cols,"%d",Size);	/* get printed width of Size */
    Sizew=strlen(Cols);
    if (Sizew>Width) Width=Sizew;
    sprintf(Entrys,"%%-%d.%df ",Width,Prec);	/* make format strs */
    sprintf(Cols,"%%-%dd ",Width);
    sprintf(Rows,"%%%dd | ",Sizew);
    Items=(Linewidth-Sizew-3)/(Width+1);	/* columns per chunk */

    /* main cycle: print chunks of the matrix */
    for (Jbeg=0,Chunk=(Size-1)/Items+1; Chunk>0; Jbeg+=Items,Chunk--)
    {
	/* underline length */
	Ulinelen=(Chunk>1)? Items*(Width+1)+Sizew+3:
		(Size-Jbeg)*(Width+1)+Sizew+3;
	/* print head line */
	for (k=0; k<Sizew+3; k++) putchar(' ');
	for (j=Jbeg; j<Size && j<Jbeg+Items; j++)
	    printf(Cols,j);	/* col numbers */
	putchar('\n');
	for (k=0; k<Ulinelen; k++) putchar('-');
	putchar('\n');

	/* print chunks for all rows */
	for (i=0; i<Size; i++)
	{
	    printf(Rows,i); 	/* row idx */
	    for (j=Jbeg; j<=i && j<Jbeg+Items; j++)
		printf(Entrys,Mat[i][j]);
	    putchar('\n');
	}

	/* chunk separator */
	putchar('\n');
	for (k=0; k<Ulinelen; k++) putchar('=');
	puts("\n");
    }
}
/* END of list_trimat */

/* list_matptr: lists the address of the matrix (i.e. the pointer Mat
 * itself) and the row pointers (the row array addresses) to 'stderr'.
 * Provided for debugging purposes.
 */
extern void list_matptr(double **Mat, int Rno)
{
    int i;
    
    fprintf(stderr, "Matrix address=%x\n", Mat); fflush(stderr);
    for (i=0; i<Rno; i++)
    {
	fprintf(stderr, "[%d]=%x\n", i, Mat[i]);
	fflush(stderr);
    }
}
/* END of list_matptr */

/* ---- SQUARE MATRICES ---- */

/* alloc_sqmat: allocates space for a square matrix (Size*Size).
 Returns the pointer to the matrix or NULL if alloc failed. */

extern Sqmat_ alloc_sqmat(int Size)
{
    Sqmat_ Mat;
    int i;
  
    Mat=(double **) calloc(Size,sizeof(double *));
    if (Mat!=NULL)
	for (i=0; i<Size; i++)
	{
	    Mat[i]=(double *)calloc(Size,sizeof(double));
	    if (Mat[i]==NULL)
	    {	/* clean up already allocated rows */
		free_matrix(Mat,i);
		return(NULL);
	    }
	}
    return(Mat);
}
/* END of alloc_sqmat */

/* list_sqmat: lists Mat to stdout with entries occupying Width chars,
 Prec digits precision. If a row takes up more than Linewidth chars,
 then the matrix is cut up nicely. */

extern void list_sqmat(Sqmat_ Mat, int Size, int Linewidth,
	int Width, int Prec)
{
    char Entrys[10],Cols[8],Rows[8];	/* format strings and len determ */
    int i,j,k,Sizew,Chunk,Jbeg,Items,Ulinelen;

    sprintf(Cols,"%d",Size);	/* get printed width of Size */
    Sizew=strlen(Cols);
    if (Sizew>Width) Width=Sizew;
    sprintf(Entrys,"%%-%d.%df ",Width,Prec);	/* make format strs */
    sprintf(Cols,"%%-%dd ",Width);
    sprintf(Rows,"%%%dd | ",Sizew);
    Items=(Linewidth-Sizew-3)/(Width+1);	/* columns per chunk */

    /* main cycle: print chunks of the matrix */
    for (Jbeg=0,Chunk=(Size-1)/Items+1; Chunk>0; Jbeg+=Items,Chunk--)
    {
	/* underline length */
	Ulinelen=(Chunk>1)? Items*(Width+1)+Sizew+3:
		(Size-Jbeg)*(Width+1)+Sizew+3;
	/* print head line */
	for (k=0; k<Sizew+3; k++) putchar(' ');
	for (j=Jbeg; j<Size && j<Jbeg+Items; j++)
	    printf(Cols,j);	/* col numbers */
	putchar('\n');
	for (k=0; k<Ulinelen; k++) putchar('-');
	putchar('\n');

	/* print chunks for all rows */
	for (i=0; i<Size; i++)
	{
	    printf(Rows,i); 	/* row idx */
	    for (j=Jbeg; j<Size && j<Jbeg+Items; j++)
		printf(Entrys,Mat[i][j]);
	    putchar('\n');
	}

	/* chunk separator */
	putchar('\n');
	for (k=0; k<Ulinelen; k++) putchar('=');
	puts("\n");
    }
}
/* END of list_sqmat */

/* ---- LU-DECOMPOSITION ---- */

/* lu_decomp: performs an LU-decomposition in place on the n*n matrix
 * A. Based on partial pivoting: the row permutations are done in Perm[]
 * (allocated within!) and will be used by lu_solve().
 * If **Perm==NULL, then the permutation vector will be used
 * internally and will be freed before return. This option is
 * used when only the determinant is calculated from the LU-decomposition.
 * Return value: the sign of the determinant of the permutation
 * matrix (+/-1) or 0 if A is singular or n<=0.
 */
int lu_decomp(Sqmat_ A, int n, int **Perm)
{
    register int i, j, k, imax, Psign=1;
    int *Idx=NULL;
    register double Large, Pivot, Tmp, Tmp2;
    double *Scal=NULL;
    
    /* check and array initialisation */
    if (n<=0)
    {
	fprintf(stderr, "? lu_decomp(): invalid matrix size %d\n", n);
	if (Perm!=NULL) *Perm=NULL; 
	return(0);
    }
    Idx=(int *) calloc(n, sizeof(int));	    /* permutation vector */
    Scal=(double *) calloc(n, sizeof(double));	/* implicit scaling array */
    
    /* get implicit scaling: if a row contains 0-s only,
     * then the matrix is singular which will be indicated
     * by setting Psign=0. Precision is controlled by
     * the constant LU_EPSILON (see Definitions above).
     */
    for (i=0; Psign && i<n; i++)
    {
	Large=0.0;
	for (j=0; j<n; j++)
	    if ((Tmp=fabs(A[i][j]))>Large) Large=Tmp;
	if (Large<LU_EPSILON)	/* (almost) singular */
	{ Psign=0; break; }
	Scal[i]=1.0/Large;
    }
    
    /* loop over columns */
    for (j=0; Psign && j<n; j++)
    {
	for (i=0; i<j; i++)
	{
	    Tmp=A[i][j];
	    for (k=0; k<i; k++) Tmp-=A[i][k]*A[k][j];
	    A[i][j]=Tmp;
	}
	
	/* find largest pivot */
	Large=0.0;
	for (i=j; i<n; i++)
	{
	    Tmp=A[i][j];
	    for (k=0; k<j; k++) Tmp-=A[i][k]*A[k][j];
	    A[i][j]=Tmp;
	    
	    if ((Tmp2=Scal[i]*fabs(Tmp))>=Large)    /* best so far */
	    { Large=Tmp2; imax=i; }
	}
	
	/* interchange rows? */
	if (j!=imax)
	{
	    for (k=0; k<n; k++)	    /* not too efficient copy */
	    {
		Tmp=A[imax][k]; A[imax][k]=A[j][k]; 
		A[j][k]=Tmp;
	    }
	    Psign*=(-1);    /* parity change */
	    Scal[imax]=Scal[j];
	}
	Idx[j]=imax;
	
	/* get the pivot */
	Pivot=A[j][j];
	if (fabs(Pivot)<LU_EPSILON)	/* singularity */
	{ Psign=0; break; }
	
	/* divide by the pivot */
	if (j<n-1)
	    for (i=j+1; i<n; i++) A[i][j]/=Pivot;
    }	/* for j */
    
    free(Scal); 
    if (Perm!=NULL) *Perm=Idx; else free(Idx);
    return(Psign);
}
/* END of lu_decomp */

/* lu_det: calculates the determinant of the n x n LU-decomposed
 * square matrix Lu. Psign is the permutation sign returned by
 * lu_decomp().
 */
double lu_det(const Sqmat_ Lu, int Psign, int n)
{
    register int i;
    register double Det, Aii;
    
    if (!Psign) return(0.0);	/* matrix is singular */
    
    Det=0.0;
    for (i=0; i<n; i++)
    {
	if ((Aii=Lu[i][i])<0.0) Psign*=(-1);   /* save sign */
	Det+=log(fabs(Aii));	/* sum logarithms to avoid overflow */
    }
    Det=Psign*exp(Det);
    return(Det);
}
/* END of lu_det */

/* lu_solve: solves the linear equation A*x=b (A is n*n, x,b are n long).
 * A is supposed to have been LU-decomposed by lu_decomp() above and
 * the row permutation is stored in Perm[]. b[] is the "right-hand-side"
 * vector which contains the solution on return.
 */
void lu_solve(const Sqmat_ A, const int Perm[], double b[], int n)
{
    register int i, j, ip;
    register double Tmp;

    /* permute forward */
    for (i=0; i<n-1; i++)
	if ((ip=Perm[i])!=i)
	{ Tmp=b[ip]; b[ip]=b[i]; b[i]=Tmp; }
    
    /* forward substitution */
    for (i=0; i<n; i++)
    {
	Tmp=b[i];
	for (j=0; j<i; j++) Tmp-=A[i][j]*b[j];
	b[i]=Tmp;
    }
    
    /* back substitution */
    for (i=n-1; i>=0; i--)
    {
	Tmp=b[i];
	for (j=i+1; j<n; j++) Tmp-=A[i][j]*b[j];
	b[i]=Tmp/A[i][i];
    }
}
/* END of lu_solve */

#undef LU_EPSILON

/* ==== END OF FUNCTIONS matrix.c ==== */

/* ==== FUNCTIONS ql.c ==== */

/* Eigenvalues and eigenvectors of real symmetric matrices by
 Housholder tridiagonalisation and QL-transformation.
 Adapted from Numerical Recipes. The "core" routines use 
 the 1..N indexing convention but the "shell" is ordinary 
 0..N-1 C-style. Standalone version is eigenql.c */

/* ANSI C, Iris Indigo IRIX 4.0.5, 20. Nov. 1992. Andris Aszodi */

/* ---- HEADER FILES ---- */

//#include "incl/ql.h"

/* ---- TYPES ---- */

/* type of function returning int: demanded by qsort() */
typedef int (*Compfnc_)(const void*,const void*);

/* type for eigenvalue sorting */
typedef struct
{
    double Eig;	/* an eigenvalue along with its */
    int Idx;	/* original position */
} Eigidx_ ;

/* ---- PROTOTYPES ---- */

/* eigen_ql(): prototype in "ql.h" */

static int eig_cmp(const Eigidx_ *E1, const Eigidx_ *E2);
static void tred2(double **a, int n, double *d, double *e);
static int tqli(double *d, double *e, int n, double **z, int Itno);

/* ---- FUNCTIONS ---- */

/* ---- SHELL ---- */

/* eigen_ql: a 'shell' function driving the Housholder and QL routines.
 Takes a lower triangular matrix Mat as input (with size Size), and
 produces the eigenvalues in Eval and the eigenvectors in Evec.
 Eval and Evec are assumed to be allocated outside with the
 correct sizes! Index shifts are performed
 to hack around the [1..N] convention of Numerical Recipes.
 Return value: 0 if OK, 1 if iteration limit was exceeded in tqli(). */ 

extern int eigen_ql(Trimat_ Mat, int Size, double *Eval, Matrix_ Evec)
{
    const int ITERNO=30;

    int i,j,k,Err;
    Sqmat_ Qmat;
    double *Diag2;
    Eigidx_ *Evs;

    /* allocate full square matrix: input and raw eigenvectors */
    Qmat=alloc_sqmat(Size);

    /* copy values from lower triangle Mat */
    for (i=0; i<Size; i++)
    {
	Qmat[i][i]=Mat[i][i];
	for (j=0; j<i; j++)
	    Qmat[i][j]=Qmat[j][i]=Mat[i][j];
    }

    /* shift addresses so that [1..N] indexing be valid */
    for (i=0; i<Size; i++) Qmat[i]--;
    Qmat--;

    /* create array for 2nd diagonal of tridiagonal matrix */
    Diag2=(double *) calloc(Size, sizeof(double));

    /* shift their address as well for [1..N] indexing */
    Eval--; Diag2--;

    /* perform Housholder tridiagonalisation */
    tred2(Qmat,Size,Eval,Diag2);

    /* apply shifted QL-transforms: get eigenvalues in Diag,
     eigenvectors in Qmat */
    Err=tqli(Eval,Diag2,Size,Qmat,ITERNO);
    if (Err) fprintf(stderr,"Iteration limit (%d) is exceeded\n",ITERNO);

    /* shift back addresses for C-style indexing */
    Eval++; Diag2++; free(Diag2);		/* not even needed */
    for (i=1; i<=Size; i++) Qmat[i]++;
    Qmat++;		/* raw eigenvectors in columns */

    /* sort eigenvalues in decreasing order */
    Evs=(Eigidx_ *) calloc(Size,sizeof(Eigidx_));
    for (i=0; i<Size; i++)
    {
	/* copy eigenvalues (rounded to 0.0 if necessary) and index */
	Evs[i].Eig=RND0(Eval[i]);
	Evs[i].Idx=i;
    }
    /* Quicksort from stdlib */
    qsort(Evs,Size,sizeof(Eigidx_),(Compfnc_)eig_cmp);
    /* permute eigenvectors according to their eigenvalues:
     and list them as rows */
    for (i=0; i<Size; i++)
    {
	Eval[i]=Evs[i].Eig;	/* copy back eigenvalues */
	k=Evs[i].Idx;
	for (j=0; j<Size; j++)
	    Evec[i][j]=Qmat[j][k];	/* copy eigenvectors */
    }

    /* output and cleanup */
    free(Evs);
    free_matrix(Qmat,Size);
    return(Err);
}
/* END of eigen_ql */

/* ---- AUXILIARY ROUTINES TO SHELL ---- */

/* eig_cmp: compares two Eigidx structs for qsort(). */

static int eig_cmp(const Eigidx_ *E1, const Eigidx_ *E2)
{
    return((E2->Eig>E1->Eig)? 1: (E2->Eig<E1->Eig)? -1: 0);
}
/* END of eig_cmp */

/* ---- EIGENROUTINE CORE ---- */

/* tred2: Housholder tridiagonalisation. a is a real, symmetric matrix
 (size n*n). The main diagonal of the tridiag. output is returned in d,
 the second diagonal in e with e[1]==0.0. On return, a contains the
 transformation matrix "Q". [1..N] indexing convention is used 
 throughout!
 Algorithm and implementation from Numerical Recipes. A.A. has after-
 edited the function head to look like ANSI, and inserted 'register'
 local vars. */

static void tred2(double **a, int n, double *d, double *e)
{
    register int l,k,j,i;
    register double scale,hh,h,g,f;

    for (i=n;i>=2;i--) {
	l=i-1;
	h=scale=0.0;
	if (l > 1) {
	    for (k=1;k<=l;k++)
		scale += fabs(a[i][k]);
	    if (scale < EPSILON)
		e[i]=a[i][l];
	    else {
		for (k=1;k<=l;k++) {
		    a[i][k] /= scale;
		    h += a[i][k]*a[i][k];
		}
		f=a[i][l];
		g = (RND0(f)>0.0) ? -sqrt(h) : sqrt(h);
		e[i]=scale*g;
		h -= f*g;
		a[i][l]=f-g;
		f=0.0;
		for (j=1;j<=l;j++) {
		/* Next statement can be omitted if eigenvectors not wanted */
		    a[j][i]=a[i][j]/h;
		    g=0.0;
		    for (k=1;k<=j;k++)
			g += a[j][k]*a[i][k];
		    for (k=j+1;k<=l;k++)
			g += a[k][j]*a[i][k];
		    e[j]=g/h;
		    f += e[j]*a[i][j];
		}
		hh=f/(h+h);
		for (j=1;j<=l;j++) {
		    f=a[i][j];
		    e[j]=g=e[j]-hh*f;
		    for (k=1;k<=j;k++)
			a[j][k] -= (f*e[k]+g*a[i][k]);
		}
	    }
	} else
	    e[i]=a[i][l];
	d[i]=h;
    }
    /* Next statement can be omitted if eigenvectors not wanted */
    d[1]=0.0;
    e[1]=0.0;
    /* Contents of this loop can be omitted if eigenvectors not
	    wanted except for statement d[i]=a[i][i]; */
    for (i=1;i<=n;i++) {
	l=i-1;
	if (RND0(d[i])!=0.0) {	/* !=0.0 added */
	    for (j=1;j<=l;j++) {
		g=0.0;
		for (k=1;k<=l;k++)
		    g += a[i][k]*a[k][j];
		for (k=1;k<=l;k++)
		    a[k][j] -= g*a[k][i];
	    }
	}
	d[i]=(fabs(a[i][i])<EPSILON)? 0.0: a[i][i];
	a[i][i]=1.0;
	for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
    }
}
/* END of tred2 */

/* tqli: QL algorithm with implicit shifts on tridiagonal matrices.
 The main diagonal is in d, the second diagonal is in e, with e[1]
 ignored. Size is n. If d and e were obtained from a general symmetric
 real matrix by Housholder transformation by tred2(), then z should
 contain the transformation matrix "Q" on input; otherwise it should
 be the unit matrix. On output, d contains the eigenvalues and z the
 eigenvectors, with the k-th column corresponding to the k-th eigenvalue.
 Itno supplies the maximum allowable no. of iterations (an addition
 by A.A.). Return value: 0 if OK, 1 if iteration limit has been
 exceeded (added by A.A. to replace the nrerror() error message function
 originally used in Numerical Recipes routines).
 [1..N] indexing convention is used throughout!
 Algorithm and implementation from Numerical Recipes. A.A. has after-
 edited the function head to look like ANSI, and inserted 'register'
 local vars. */

#define SIGN(a,b) ()

static int tqli(double *d, double *e, int n, double **z, int Itno)
{
    register int m,l,iter,i,k;
    register double s,r,ra,p,g,f,c,b;
    register float dd;

    for (i=2;i<=n;i++) e[i-1]=e[i];
    e[n]=0.0;
    for (l=1;l<=n;l++) {
	iter=0;
	do {
	    for (m=l;m<=n-1;m++) {
		dd=(float)(fabs(d[m])+fabs(d[m+1]));
		if ((float)fabs(e[m])+dd ==dd) break; 
	    }
	    if (m != l) {
		if (iter++ >= Itno)	/* too many iters */
		    return(1);
		g=(d[l+1]-d[l])/(2.0*e[l]);
		r=sqrt((g*g)+1.0);
		ra=(RND0(g)<0.0) ? -fabs(r) : fabs(r);
		g=d[m]-d[l]+e[l]/(g+ra);
		s=c=1.0;
		p=0.0;
		for (i=m-1;i>=l;i--) {
		    f=s*e[i];
		    b=c*e[i];
		    if (fabs(f) >= fabs(g)) {
			c=g/f;
			r=sqrt((c*c)+1.0);
			e[i+1]=f*r;
			c *= (s=1.0/r);
		    } else {
			s=f/g;
			r=sqrt((s*s)+1.0);
			e[i+1]=g*r;
			s *= (c=1.0/r);
		    }
		    g=d[i+1]-p;
		    r=(d[i]-g)*s+2.0*c*b;
		    p=s*r;
		    d[i+1]=g+p;
		    g=c*r-b;
		    /* Next loop can be omitted if eigenvectors not wanted */
		    for (k=1;k<=n;k++) {
			f=z[k][i+1];
			z[k][i+1]=s*z[k][i]+c*f;
			z[k][i]=c*z[k][i]-s*f;
		    }
		}
		d[l]=d[l]-p;
		e[l]=g;
		e[m]=0.0;
	    }
	} while (m != l);
    }
    return(0);
}
/* END of tqli */

/* ==== END OF FUNCTIONS ql.c ==== */
