#include "util.hpp"

double sign ( double x ) { if (x<0) return -1.0; else return 1.0; }
float  sign (  float x ) { if (x<0) return -1.0; else return 1.0; }
int    sign (    int x ) { if (x<0) return -1;   else return 1;   }

// random float (0..1)
float randf () { return (float)(rand()%100000)*0.00001; }

// io functions

int read_line( FILE *file, char *string )
{	char	c;
	int	i=0;
	*string = 0;
	while(c=getc(file)) {
		/* printf("%d >%c<\n", c,c); */
		if (feof(file)) return -i-1;
		if (c=='\n') return i;
		string[i] = c;
		string[++i] = (char)0;
	}
}
int next_line ( FILE *file )
{	char	c;
	while(c=getc(file)) {
		if (feof(file)) return 0;
		if (c=='\n') return 1;
	}
}

// timing functions

/*
timespec start,finish;
	clock_gettime(CLOCK_REALTIME, &start);
	stuff
	clock_gettime(CLOCK_REALTIME, &finish);
	runtime = timedifL(start,finish);
			 F or D for seconds
*/

#define NANO  1000000000
#define NANOf 1000000000.0

long int timedif (struct timespec start, struct timespec end)
{ // not used outside here
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = NANO+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return NANO*temp.tv_sec + temp.tv_nsec;
}

long int timedifL (struct timespec start, struct timespec end) {
	return timedif(start,end);
}

float timedifF (struct timespec start, struct timespec end) {
	return (float)timedif(start,end)/NANOf;
}

double timedifD (struct timespec start, struct timespec end) {
	return (double)timedif(start,end)/NANOf;
}

void timeout ( long int nanosecs ) {
struct  timespec wait;
int	seconds;
	if ( nanosecs > NANO ) {
		seconds = (int)(nanosecs/NANO);
		nanosecs -= seconds*NANO;
	} else {
		wait.tv_sec = 0;
	}
	wait.tv_nsec = nanosecs;
	nanosleep(&wait,NULL);
}

void timeout ( int seconds, long int nanosecs ) {
struct  timespec wait;
	wait.tv_sec = seconds;
	wait.tv_nsec = nanosecs;
	nanosleep(&wait,NULL);
}

// sorting functions

int sort4min ( float a, float b, float c, float d ) {
	if ( a<b && a<c && a<d ) return 1;
	if ( b<a && b<c && b<d ) return 2;
	if ( c<a && c<b && c<d ) return 3;
	if ( d<a && d<b && d<c ) return 4;
	return 0;
}

int sort4max ( float a, float b, float c, float d ) {
	if ( a>b && a>c && a>d ) return 1;
	if ( b>a && b>c && b>d ) return 2;
	if ( c>a && c>b && c>d ) return 3;
	if ( d>a && d>b && d>c ) return 4;
	return 0;
}

int sort4min ( int a, int b, int c, int d ) {
	if ( a<b && a<c && a<d ) return 1;
	if ( b<a && b<c && b<d ) return 2;
	if ( c<a && c<b && c<d ) return 3;
	if ( d<a && d<b && d<c ) return 4;
	return 0;
}

int sort4max ( int a, int b, int c, int d ) {
	if ( a>b && a>c && a>d ) return 1;
	if ( b>a && b>c && b>d ) return 2;
	if ( c>a && c>b && c>d ) return 3;
	if ( d>a && d>b && d>c ) return 4;
	return 0;
}

/* SHELL SORTS ON (p/f/i)a WITH POINTERS p (a IS UNCHANGED), BIGGEST TO TOP */

void sort ( Pairs *pa, float *fa, int *ia, int *p, int n, bool init_pointers )
{
int	i, j, m;
int	reverse = 0; // -ve n reverses sorting order
int	pair = 0, floating = 0, integer = 0;

	if ( n < 0 ) { n = -n; reverse = 1; }
	if ( n == 0 ) return;					/* DEAL WITH TRIVIAL CASES */
	if ( n == 1 ) {
		p[0] = 0;
		return;
	}
	if (ia) integer = 1;					/* SET NUMBER MODE */
	if (fa) floating = 1;
	if (pa) pair = 1;
	if( init_pointers ) for (i=0; i<n; i++) p[i] = i;	/* INITIALISE THE POINTERS */
	m = n;
	while (1) { 						/* LOOP OF DECREASING SWAP SPAN	*/
		m = m/2;
		if(m == 0) return;
		for ( i = 0; i < n-m; i++ ) {			/* LOOP OVER SWAPS */
			j = i;					/* LOOP TO BUBBLE UP */
			while (1) { int	jp;
				if (reverse) {
					if ( pair   && pa[p[j+m]].s >= pa[p[j]].s ) break;
					if ( floating && fa[p[j+m]] >= fa[p[j]] ) break;
 					if ( integer  && ia[p[j+m]] >= ia[p[j]] ) break;
				} else {
					if ( pair   && pa[p[j+m]].s <= pa[p[j]].s ) break;
					if ( floating && fa[p[j+m]] <= fa[p[j]] ) break;
 					if ( integer  && ia[p[j+m]] <= ia[p[j]] ) break;
				}
				jp = p[j];
				p[j] = p[j+m];			/* EXCHANGE POINTERS */
				p[j+m] = jp;
				j = j-m;			/* RETURN TO TOP OF LIST */
				if( j < 0 ) break;		/* ONLY IF ROOM */
			}
		}
	}
}
void sort ( Pairs *a, int *p, int n ) { sort(a,0,0,p,n,1); }
void sort ( float *a, int *p, int n ) { sort(0,a,0,p,n,1); }
void sort ( int   *a, int *p, int n ) { sort(0,0,a,p,n,1); }

// SHELL SORT with reordering

void sort ( float *a, int n )
{
int	i, j, m;
	if ( n <= 1 ) return;					/* DEAL WITH TRIVIAL CASES */
	m = n;
	while (1) { 						/* LOOP OF DECREASING SWAP SPAN	*/
		m = m/2;
		if(m == 0) return;
		for ( i = 0; i < n-m; i++ ) {			/* LOOP OVER SWAPS */
			j = i;					/* LOOP TO BUBBLE UP */
			while (1) { float b;
 				if ( a[j+m] >= a[j] ) break;
				b = a[j];
				a[j] = a[j+m];			/* EXCHANGE VALUES */
				a[j+m] = b;
				j = j-m;			/* RETURN TO TOP OF LIST */
				if( j < 0 ) break;		/* ONLY IF ROOM */
			}
		}
	}
}

void sort ( int *a, int n )
{
int	i, j, m;
	if ( n <= 1 ) return;					/* DEAL WITH TRIVIAL CASES */
	m = n;
	while (1) { 						/* LOOP OF DECREASING SWAP SPAN	*/
		m = m/2;
		if(m == 0) return;
		for ( i = 0; i < n-m; i++ ) {			/* LOOP OVER SWAPS */
			j = i;					/* LOOP TO BUBBLE UP */
			while (1) { int b;
 				if ( a[j+m] >= a[j] ) break;
				b = a[j];
				a[j] = a[j+m];			/* EXCHANGE VALUES */
				a[j+m] = b;
				j = j-m;			/* RETURN TO TOP OF LIST */
				if( j < 0 ) break;		/* ONLY IF ROOM */
			}
		}
	}
}

void sort ( short *a, int n )
{
int	i, j, m;
	if ( n <= 1 ) return;					/* DEAL WITH TRIVIAL CASES */
	m = n;
	while (1) { 						/* LOOP OF DECREASING SWAP SPAN	*/
		m = m/2;
		if(m == 0) return;
		for ( i = 0; i < n-m; i++ ) {			/* LOOP OVER SWAPS */
			j = i;					/* LOOP TO BUBBLE UP */
			while (1) { int b;
 				if ( a[j+m] >= a[j] ) break;
				b = a[j];
				a[j] = a[j+m];			/* EXCHANGE VALUES */
				a[j+m] = b;
				j = j-m;			/* RETURN TO TOP OF LIST */
				if( j < 0 ) break;		/* ONLY IF ROOM */
			}
		}
	}
}

// bipartite graph matching

void makepair ( int, int, int**, int**, Pairs* );

float pairup ( float **mat, int n0, int n1, Pairs *pairs, int *n)
// takes the score matrix <mat> (n0xn1, high=good) and returns <n> matched pairs
{
int i,j;
float   *ranks, score;
int     **rank0, **rank1;
        ranks = new float[max(n0,n1)];
        rank0 = new int*[n0];
        for (i=0; i<n0; i++) rank0[i] = new int[n1];
        rank1 = new int*[n1];
        for (i=0; i<n1; i++) rank1[i] = new int[n0];
        for (i=0; i<n0; i++) {
                for (j=0; j<n1; j++) ranks[j] = mat[i][j]+randf()*0.001;
                sort(0,ranks,0,rank0[i],n1,1);
        }
        for (i=0; i<n1; i++) {
                for (j=0; j<n0; j++) ranks[j] = mat[j][i]+randf();0.001;
                sort(ranks,rank1[i],n0);
        }
/*
        for (i=0; i<n0; i++) { for (j=0; j<n1; j++) printf("%5.1f", mat[i][j]); NL } NL
        for (i=0; i<n0; i++) { for (j=0; j<n1; j++) printf("%3d", rank0[i][j]); NL } NL
        for (i=0; i<n1; i++) { for (j=0; j<n0; j++) printf("%3d", rank1[i][j]); NL } NL
*/
        score = 0.0;
        *n = min(n0,n1);
        if (n0<n1) {
                makepair(n0,n1,rank0,rank1,pairs);
                for (i=0; i<(*n); i++) { int m;
                        pairs[i].s = mat[pairs[i].a][pairs[i].b];
			m = pairs[i].a; pairs[i].a=pairs[i].b; pairs[i].b=m;
                        score += pairs[i].s;
                }
        } else {
                makepair(n1,n0,rank1,rank0,pairs);
                for (i=0; i<(*n); i++) {
                        pairs[i].s = mat[pairs[i].b][pairs[i].a];
                        score += pairs[i].s;
                }
        }
        return score;
}

void makepair ( int n0, int n1, int **rank0, int **rank1, Pairs *pairs )
// bipartite graph matching algorithm (not used outside here)
{
int     i,j,k,n,*at = new int[n0];
        for (i=0; i<n0; i++) { at[i] = 0; pairs[i].a = i; pairs[i].b = -1; }
Pi(n0) Pi(n1) NL
        n = 0;
        while (n<n0) {
                for (i=0; i<n0; i++) {
                        if (pairs[i].b >= 0) continue;
                        pairs[i].b = rank0[i][at[i]];
                        for (j=0; j<n0; j++) {
                                if (j==i) continue;
                                if (pairs[i].b == pairs[j].b) { int b = pairs[i].b;
                                /* partner already allocated */
                                        for (k=0; k<n1; k++) {
                                                if (rank1[b][k]==pairs[i].a) {
                                                /* i gets to keep it */
                                                        pairs[j].b = -1;
                                                        n--;
                                                        break;
                                                }
                                                if (rank1[b][k]==pairs[j].a) {
                                                /* j gets to keep it */
                                                        pairs[i].b = -1;
                                                        n--;
                                                        break;
                                                }
                                        }
                                        break;
                                }
                        }
                        at[i]++; n++;
                }
        }
}
