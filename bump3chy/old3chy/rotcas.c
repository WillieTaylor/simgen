/*
cc -O rotcas.c -o rotcas util/wt/util.o util/wt/geom.o util/wt/sort.o -lm -m32

*/
#include "util/wt/incl/util.h"
#include "util/wt/incl/geom.h"

#define NALLOC 1500
#define BUMP 0.5
#define CYCLES 500

Vec    res[NALLOC];
Vec    ret[NALLOC];
Vec    ray[NALLOC];

float	bond;
float	packs();

main(argc,argv)
int argc; char *argv[];
{
int	i, j, n, len, bumps;
char	line[225];
FILE	*pdb, *out;
int	beg, end, cycles, check;
float	rx, ry, rz;
Vec    	x,y,z, cog, zero;
	vinit(&zero);
	vinit(&x); x.x = 1.0;
	vinit(&y); y.y = 1.0;
	vinit(&z); z.z = 1.0;
	pdb = fopen(argv[1],"r");
	sscanf(argv[2],"%f", &rx);
	sscanf(argv[3],"%f", &ry);
	sscanf(argv[4],"%f", &rz);
	rx = rx*PI/180.0;
	ry = ry*PI/180.0;
	rz = rz*PI/180.0;
	len = getca(ray,pdb);
	vinit(&cog);
	for (i=1; i<=len; i++) vsum(ray[i],&cog);
	vdiv(&cog,(float)len);
	for (i=1; i<=len; i++) vsub(ray[i],cog,res+i);
	for (i=1; i<=len; i++) {
		rotate(zero,x,res+i,rx);
		rotate(zero,y,res+i,ry);
		rotate(zero,z,res+i,rz);
	}
	out = fopen("rot.out","w");
	putpdb (res,out,len,1.0);
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
 
putpdb (res,out,len,sca)
Vec    *res;
FILE    *out;
int     len;
float	sca;
{       int     i = 0, n = 0;
        for (i=1; i<=len; i++) {
                	fprintf(out,"ATOM%7d  CA  GLY%6d     %7.3f %7.3f %7.3f   0.0   0.0\n",
                        i, i, sca*res[i].x, sca*res[i].y, sca*res[i].z);
        }
        fprintf(out,"TER\n");
}
