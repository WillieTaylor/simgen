#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <math.h>
#include <time.h>
#include <string.h>
#include <pthread.h>

using namespace std;

#define Pc(c)   cout<<" "#c" = "<<c;
#define Ps(s)   cout<<" "#s" = "<<s;
#define Pi(i)   cout<<" "#i" = "<<i;
#define Pr(r)   cout<<" "#r" = "<<r;
#define	Pv(v)	cout<<" "#v" = "<<v.x<<" "<<v.y<<" "<<v.z;
#define	Pl(s)	cout<<" "#s" = "<<s.A.x<<" "<<s.A.y<<" "<<s.A.z<<" --- "<<s.B.x<<" "<<s.B.y<<" "<<s.B.z;
#define	Pm(m)	cout<<" "#m" = "<<m.A.x<<" "<<m.A.y<<" "<<m.A.z<<", "<<m.B.x<<" "<<m.B.y<<" "<<m.B.z<<", "<<m.C.x<<" "<<m.C.y<<" "<<m.C.z;
#define NL	cout<<"\n";
#define PT	cout<<"\t";
#define Px	cout<<"here\n";
#define Pt(t)	cout<<" "#t" ";
#define DO(i,n)	 for (int i=0; i<n; i++)
#define DO1(i,n) for (int i=1; i<=n; i++)
#define LOOP	 while (1)
#define TRUE 1
#define FALSE 0

#define TEST(obj)    if (obj == NULL) { printf("malloc fail for " #obj "\n"); exit(1); }

const float PI = 3.14159265359f;
const float twoPI = PI*2.0f;
const float halfPI = PI*0.5f;
const float NOISE = 0.00001f;
const float BIG = 10000000.0f;

// get sign
double sign ( double );
float  sign ( float  );
int    sign ( int    );

// random float (0..1)
float randf ();

// io functions
int read_line( FILE*, char* );
int next_line ( FILE* );

int acid2num ( char* );

// timing functions

long int timedifL ( struct timespec, struct timespec );
float	 timedifF ( struct timespec, struct timespec );
double	 timedifD ( struct timespec, struct timespec );

void timeout ( long );
void timeout ( int, long );

// sorting functions

int sort4min ( float, float, float, float );
int sort4max ( float, float, float, float );
int sort4min ( int, int, int, int );
int sort4max ( int, int, int, int );

typedef struct { int a,b; float s; char c; } Pairs;

// shell sorts on (p/f/i)a with pointers p (a is unchanged), biggest to top
void sort ( Pairs*, float*, int*, int*, int, bool);
void sort ( Pairs*, int*, int );
void sort ( float*, int*, int );
void sort ( int  *, int*, int );

// shell sort (changes a)
void sort ( float*, int );
void sort (   int*, int );
void sort ( short*, int );

// bipartite graph matching
float pairup ( float**, int, int, Pairs*, int );
