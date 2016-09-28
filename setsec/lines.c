#include "util/wt/incl/util.h"
#include "util/wt/incl/geom.h"
#include "util/aa/incl/matrix.h"

#define BREAK 10.0

float	linefit();
float	recp_area();
float	recp_ends();
float	recp_end();
double	supermac();

float trieq_bal(Trimat_ Metric, int Rno, float *Diagshf);
int metric_project(Trimat_ Metric, int Rno, Sqmat_ Xyz,
        float Diagshf, float Equ, int Oldim, float *Projqual);

lines (gap,beta,a,c,b,len,ends,dens,vars)
/* takes alpha-carbon coordinates in a[] (1...N) and
   returns the line end-points in b[] (with termini in ends[])
   native coordinates (not smoothed) are in c[] */
float gap, beta;
Vec *a,*b,*c; int len, *ends; float *dens, *vars;
{
int	i, j, n, cycles = 2;
int	m, w, topm, topw;
float	tops;
float	**mat, **sum, *emax;
int	**top[2], *edge[2];
	mat = (float**)malloc(sizeof(float*)*(len+1)); TEST(mat)
	sum = (float**)malloc(sizeof(float*)*(len+1)); TEST(sum)
	top[0] = (int**)malloc(sizeof(int*)*(len+1));  TEST(top[0])
	top[1] = (int**)malloc(sizeof(int*)*(len+1));  TEST(top[1])
	for (i=0; i<=len; i++) {
		mat[i] = (float*)malloc(sizeof(float)*(len+1)); TEST(mat[i])
		sum[i] = (float*)malloc(sizeof(float)*(len+1)); TEST(sum[i])
		top[0][i] = (int*)malloc(sizeof(int)*(len+1));  TEST(top[0][i])
		top[1][i] = (int*)malloc(sizeof(int)*(len+1));  TEST(top[1][i])
	}
	emax = (float*)malloc(sizeof(float)*(len+1));
	edge[0] = (int*)malloc(sizeof(int)*(len+1));
	edge[1] = (int*)malloc(sizeof(int)*(len+1));
	for (m=0; m<len; m++) for (w=0; w<len; w++) mat[m][w] = sum[m][w] = 0.0;
	for (m=2; m<len-2; m++) { int skip; float d, e;
	/*     2 .. len-2 allows w=2 at each end  */
		for (w=2; w<len/2; w++) { float s;
			if (m-w < 0) break;
			if (m+w > len) break;
			mat[m][w] = -gap*(float)(w+1); /* -gap*sqrt((float)w); */
			skip = 0;
			sum[m][w] = 1.0; /* flag for OK = 1, forbidden = -1 */
			for (j=m-w; j<=m+w; j++) if (ends[j+1]) { skip = 1; sum[m][w] = -1.0; break; }
			if (skip) continue;
			mat[m][w] += linefit(a,a,m,w,0,&d,&e);
		}
	}
	for (m=0; m<len; m++) { /* NB ends and dens run 1...N */
		if (ends[m+1]) printf("%3d * ", m); else printf("%3d | ", m);
		sum[m][0] = mat[m][0] = 0.0;
		sum[m][1] = mat[m][1] = 0.0;
		if (dens[m+1] > 0.001) sum[m][1] = mat[m][1] = beta*log(1.0+dens[m+1]);
/*
		if (dens[m+1] > 0.001) sum[m][1] = mat[m][1] = beta*sqrt(dens[m+1]);
*/
		for (w=1; w<20; w++) printf("%3d", (int)(mat[m][w])); NL
	}
	for (w=2; w<len/2; w++) {
		for (m=2; m<len-2; m++) {
			if (m-w < 0) continue;
			if (m+w > len) continue;
			if (w==2 && sum[m][1] > 0.0) {
					sum[m][w] = mat[m][w]+sum[m][1]+sum[m-1][1]+sum[m+1][1];
			} else {
					sum[m][w] = sum[m-1][w-1]+sum[m+1][w-1]-sum[m][w-2]+mat[m][w]+mat[m][w-1];
			}
			for (j=m-w; j<=m+w; j++) { if (ends[j+1]) sum[m][w] = -1.0; }
		}
	}
	NLL
	for (m=0; m<len; m++) {
		if (ends[m+1]) printf("%3d * ", m); else printf("%3d | ", m);
		for (w=1; w<20; w++) {
			if (sum[m][w] > 0.0) printf("%4d", (int)sum[m][w]);
					else printf("    ");
		} NL
	}
	for (m=0; m<len; m++) {
		emax[m] = 0.0;
		edge[0][m] = edge[1][m] = 0;
		for (w=0; w<len; w++) top[0][m][w] = top[1][m][w] = -1;
	}
	tops = 0.0;
	for (n=3; n<len; n++) {
		m = n;
		for (w=1; w<len/2; w++)
		{ float max;
		  int	p, q, r;
			if (m-w < 0) break;
			if (m+w > len) break;
			max = 0.0;
			r = i = m-w-2; /* -2 forces one gap between segments */
			r = i = m-w-4; /* -3 forces 2 gaps  between segments */
			p = q = j = 0;
			while (i>0 && j<len/2) {
				/* scan leading edge for max */
				if (sum[i][j] > max) {
					if ( j==1 && mat[i][j]<1.0) {
						/* exclude non B segs of length 3 */
					} else {
						max = sum[i][j];
						p = i; q = j;
					}
				}
				i--; j++;
			}
			emax[r] = max;
			edge[0][r] = p;
			edge[1][r] = q;
			max += sum[m][w];
			for (i=1; i<=r; i++)
			{ float s = sum[m][w]+emax[i];
				/* scan earlier edges for better max */
				if (s > max) {
					max = s;
					p = edge[0][i];
					q = edge[1][i];
				}
			}
			top[0][m][w] = p;
			top[1][m][w] = q;
			sum[m][w] = max;
			if (sum[m][w] > tops) {
				tops = sum[m][w];
				topm = m; topw = w;
			}
			m--;
		}
	}
	n = 0;
	NLL Pr(tops) NLL
	i = topm; j = topw;
	while (i>0 && j>0) { int ii, jj;
		Pi(i) Pi(j) NL
		edge[0][n] = i;
		edge[1][n] = j;
		ii = top[0][i][j]; jj = top[1][i][j];
		i = ii; j = jj;
		n++;
	}
	NL
	m = 0;
	for (i=n-1; i>-1; i--)
	{ int	nn = edge[0][i]-edge[1][i]+1,
		cc = edge[0][i]+edge[1][i]+1;
	  float d,e;
		m++;
		linefit(a,a,edge[0][i],edge[1][i],b+m,&d,&e);
		ends[m] = nn;
		dens[m] = d;
		vars[m] = e;
		m++;
		ends[m] = cc;
		dens[m] = d;
		vars[m] = e;
		
	}
	for (i=0; i<len; i++) {
		free(mat[i]); free(sum[i]); free(top[0][i]); free(top[1][i]);
	}
/*
	free(mat); free(sum); free(top[0]); free(top[1]);
	free(emax); free(edge[0]); free(edge[1]);
*/
	return m/2;
}

smooth (a,len,cycles) Vec *a; int len, cycles;
{
int	i, j;
	Pi(cycles) Pi(len) NL
	for (i=0; i<cycles; i++) { Vec p1, p2, p3, q;
		vcopy(a[1],&p1); vcopy(a[2],&p2); vcopy(a[3],&p3);
		for (j=2; j<len; j++) {
			if (vdif(a[j+1],a[j+0]) > BREAK) continue;
			if (vdif(a[j+1],a[j+2]) > BREAK) continue;
			vinit(&q);
			vsum(p1,&q); vsum(p2,&q); vsum(p3,&q);
			vdiv(&q,3.0);
			vcopy(a[j],&p1); vcopy(a[j+1],&p2); vcopy(a[j+2],&p3);
			vcopy(q,a+j);
		}
	}
}

packout (a,b,c,d,len,m,ends,dens,leng,angs,dmat,area,pack,snet,net,id)
Vec *a,*b,*c,*d; int id,len,m, *ends;
float *dens, *leng, **angs, **dmat, **area, **pack, **net;
{
double  *ww, *ac, *bc, **va, **vb;
Sqmat_  rot1, rot2;
Mat     mat;
Vec     acnt, bcnt;
float	rms, rms1, rms2;
FILE	*out;
char	file[22];
int	i, j, n, lem = m*2;
float	*dat, *temp, *lang;
	dat = (float*)alloca(sizeof(float)*(lem+1));
	temp = (float*)alloca(sizeof(float)*(lem+1));
	lang = (float*)alloca(sizeof(float)*(lem+1));
	sprintf(file,"lines%d.pdb", id);
	out = fopen(file,"w");
	pdbput(a+1,out,len,'A',0,0);
	pdbput(b+1,out,lem,'B',0,0);
        for (i=0; i<m; i++) dat[i*2+1] = dens[i];
Pi(m) NL
	packs(b,c,ends,m,dmat,area,pack,angs,snet,net,dat,id);
        for (i=0; i<m; i++) { lang[i] = dat[i*2+2]; temp[i] = 4.0-dat[i*2+1]; }
	pdbput(c,out,m,'C',lang,temp);
        for (i=0; i<m; i++) vcopy(c[i], d+i);
	projstr(dmat,area,d,dat,m);
        rot1 = alloc_sqmat(3);
        rot2 = alloc_sqmat(3);
        ac = (double*)alloca(sizeof(double)*3);
        bc = (double*)alloca(sizeof(double)*3);
        ww = (double*)alloca(sizeof(double)*m);
        va = (double**)alloca(sizeof(double*)*m);
        vb = (double**)alloca(sizeof(double*)*m);
        for (i=0; i<m; i++) {
		ww[i] = 1.0;
                va[i] = (double*)alloca(sizeof(double)*3);
                va[i][0] = c[i].x;
                va[i][1] = c[i].y;
                va[i][2] = c[i].z;
                vb[i] = (double*)alloca(sizeof(double)*3);
                vb[i][0] = d[i].x;
                vb[i][1] = d[i].y;
                vb[i][2] = d[i].z;
	}
        rms1 = supermac(ww,va,vb,m,ac,bc,rot1);
        acnt.x = ac[0]; acnt.y = ac[1]; acnt.z = ac[2];
        bcnt.x = bc[0]; bcnt.y = bc[1]; bcnt.z = bc[2];
        printf("RMSd = %7.3f over %d atoms\n", rms1, m);
        for (i=0; i<m; i++) vb[i][2] = -d[i].z;
        rms2 = supermac(ww,va,vb,m,ac,bc,rot2);
        printf("RMSd = %7.3f over %d atoms\n", rms2, m);
	if (rms2<rms1) {
		for (i=0; i<m; i++) d[i].z *= -1.0;
		printf("Z-flip\n");
        	mat.A.x = rot2[0][0]; mat.A.y = rot2[0][1]; mat.A.z = rot2[0][2];
        	mat.B.x = rot2[1][0]; mat.B.y = rot2[1][1]; mat.B.z = rot2[1][2];
        	mat.C.x = rot2[2][0]; mat.C.y = rot2[2][1]; mat.C.z = rot2[2][2];
		rms = rms2;
	} else {
        	mat.A.x = rot1[0][0]; mat.A.y = rot1[0][1]; mat.A.z = rot1[0][2];
        	mat.B.x = rot1[1][0]; mat.B.y = rot1[1][1]; mat.B.z = rot1[1][2];
        	mat.C.x = rot1[2][0]; mat.C.y = rot1[2][1]; mat.C.z = rot1[2][2];
		rms = rms1;
	}
	if (rms > 0.0) {
        	for (i=0; i<m; i++) {
                	vsub(d[i],bcnt,d+i);
                	MmulV(&mat,d[i],d+i);
                	vadd(d[i],acnt,d+i);
		}
        }
	pdbput(d,out,m,'D',lang,temp);
	fclose(out);
	for (i=0; i<m; i++) {
		for (j=0; j<m; j++) dmat[i][j] = vdif(d[i],d[j]);
		if (dens[i] < 0.0) dens[i] = dat[i*2+1];
		leng[i] = dat[i*2+2];
Pi(i) Pr(dens[i]) Pr(leng[i]) NL
	}
}

#define CUT 30.0

packs (end,mid,res,n,dmat,smat,pack,angs,snet,net,dat,id)
Vec *end, *mid; int *res, n;
float **dmat, **smat, **pack, **angs, **snet, **net, *dat;
int id;
{
int	i, j, k, *got, nn = n*2;
float	over;
Vec	box[4];
FILE	*out[5];
	for (i=0; i<5; i++) { char name[12];
		sprintf(name,"data%c%d.dat", 'A'+id, i);
		out[i] = fopen(name,"w");
	}
Pi(n) NL
	for (i=0; i<=n; i++) for (j=0; j<=n; j++)
		{ dmat[i][j] = CUT; pack[i][j] = angs[i][j] = smat[i][j] = snet[i][j] = 0.0; }
Pi(nn) NL
for (i=1; i<nn; i+=2) { Pi(i) Pv(end[i+1]) NL }
	for (i=1; i<nn; i+=2)
	{ float	d = vdif(end[i],end[i+1]);
		dat[i+1] = d;
		if (dat[i]<0.0) dat[i] = (d+1.0)/(float)(res[i+1]-res[i]+1);
		if (dat[i] > 2.5) dat[i] = 2.5;
		Pi(i) Pi(res[i]) Pi(res[i+1]) Pr(dat[i]) Pr(dat[i+1]) NL
		vave(end[i],end[i+1],&mid[i/2]);
	}
	got = (int*)alloca(sizeof(int)*n);
	for (i=0; i<n; i++) got[i] = 1;
	for (i=1; i<nn-1; i+=2) {
		for (j=i+2; j<nn; j+=2)
		{ Vec	a, b, c, d, ab, cd, tmp;
		  float type, area, norm, d1, d2, dd, dot, dist = 0.0;
		  int	class, p = i/2, q = j/2;
Pi(i) Pi(j) NL
			vcopy(end[i],&a); vcopy(end[i+1],&b);
			vcopy(end[j],&c); vcopy(end[j+1],&d);
			vsub(a,b,&ab); vsub(c,d,&cd);
			dot = vdot(ab,cd)/(vdif(a,b)*vdif(c,d));
			dot = fabs(dot);
			if (dot>1.0) dot = 1.0;
			angs[p][q] = angs[q][p] = acos(dot);
			over = lineOline(a,b,c,d,box);
Pr(over) NL
			dd = 0.25*(vdif(a,c)+vdif(a,d)+vdif(b,c)+vdif(b,d));
			if (over<0.1 || dd>CUT ) {
				dmat[p][q] = dmat[q][p] = dd;
				smat[p][q] = smat[q][p] = 0.0;
				continue;
			}
			vave(box[0],box[1],&ab);
			vave(box[2],box[3],&cd);
			dist = vdif(ab,cd);
Pv(a) NL Pv(b) NL Pv(c) NL Pv(d) NL
			area = recp_area(box);
Pv(a) NL Pv(b) NL Pv(c) NL Pv(d) NL
Pr(area) NL
			area = recp_ends(a,b,c,d,box);
Pr(area) NL
			type = dat[i]+dat[j];
			norm = sqrt(area)/type;
			pack[p][q] = pack[q][p] = norm;
			if (norm > 0.5) {
				vave(box[0],box[1],&tmp);
				vsum(tmp,&mid[i/2]); got[i/2]++;
				vave(box[2],box[3],&tmp);
				vsum(tmp,&mid[j/2]); got[j/2]++;
			}
			dd = dist;
			if ((dat[i] > 2.0) && (dat[j] > 2.0)) { int ii,jj; float ss;
				if (net) {
					ss = 0.0;
					for (ii=res[i]; ii<=res[i+1]; ii++) 
					for (jj=res[j]; jj<=res[j+1]; jj++) ss += net[ii][jj];
				} else { ss = 10.0; }
				if ((ss>0.5) && (dd<7.5)) {
					snet[p][q] = snet[q][p] = ss; dd = 4.5;
					printf("strands %d and %d linked (bonds=%f, dist=%f ---> %f)\n",
						p, q, ss, dist, dd);
				}
			} else {
				if (dd < 10.0) dd = 10.0;
			}
			dmat[p][q] = dmat[q][p] = dd;
			smat[p][q] = smat[q][p] = area+snet[p][q];
			if (dist>CUT) continue;
			fprintf(out[0],"%d %d %d-%d %d-%d %f %f %f %f\n",
				i,j,res[i],res[i+1],res[j],res[j+1],type,dist,area,norm);
			class = (int)type;
			if (class>4) class = 4;
			if (class<1) class = 1;
			fprintf(out[class],"%d %d %d-%d %d-%d %f %f %f %f\n",
				i,j,res[i],res[i+1],res[j],res[j+1],type,dist,area,norm);
		}
	}
	for (i=0; i<n; i++) vdiv(mid+i,(float)got[i]);
	for (i=0; i<n; i++) { Pi(i) Pi(got[i]) Pv(mid[i]) NL }
}

float recp_ends(a,b,c,d,box) Vec a,b,c,d, *box;
{
int	i, ia,ib,ic,id;
float	da,db,dc,dd, area=0.0;
	da=db=dc=dd=999.9;
Pv(a) NL Pv(b) NL Pv(c) NL Pv(d) NL
	for (i=0; i<4; i++) { float p;
Pi(i) Pv(box[i]) NL
		p = vdif(box[i],a); if (p<da) { da=p; ia=i; }
		p = vdif(box[i],b); if (p<db) { db=p; ib=i; }
		p = vdif(box[i],c); if (p<dc) { dc=p; ic=i; }
		p = vdif(box[i],d); if (p<dd) { dd=p; id=i; }
Pi(i) Pi(id) Pr(p) Pr(dd) NL  
	}
Pi(ia) Pr(da) NL
Pi(ib) Pr(db) NL
Pi(ic) Pr(dc) NL
Pi(id) Pr(dd) NL
	if (da>0.0001) area += recp_end(a,box[ia],box[(ia+2)%4]);
	if (db>0.0001) area += recp_end(b,box[ib],box[(ib+2)%4]);
	if (dc>0.0001) area += recp_end(c,box[ic],box[(ic+2)%4]);
	if (dd>0.0001) area += recp_end(d,box[id],box[(id+2)%4]);
	return area;
}

projstr (mat,sat,coord,dat,len)
float   **mat, **sat;
Vec     *coord;
float	*dat;
int     len;
{
Trimat_ D, Dist, Strict;
Sqmat_ Xyz;
int     i, j, Itno;
float	smax, smin;
        D = alloc_trimat(len+1);
        Dist = alloc_trimat(len+1);
        Strict = alloc_trimat(len+1);
        Xyz = alloc_sqmat(len+1);
	smax = 0.0; smin = 999.0;
        for (i=0; i<len-1; i++) {
                for (j=i+1; j<len; j++) {
			if (sat[i][j] > smax) smax = sat[i][j];
			if (sat[i][j] < smin) smin = sat[i][j];
		}
	}
	Pr(smax) Pr(smin) NL
	smax = smax-smin;
        for (i=0; i<len; i++) for (j=0; j<len; j++) sat[i][j] = (sat[i][j]-smin)/smax;
        for (i=0; i<len; i++) { mat[i][i] = 0.0; sat[i][i] = 1.0; }
        for (i=0; i<len; i++) {
                for (j=0; j<i; j++) {
                        Strict[i][j] = sat[i][j]*2.0;
                        if (Strict[i][j] > 1.0) Strict[i][j] = 1.0;
                        if (Strict[i][j] < 0.0) Strict[i][j] = 0.0;
                        Dist[i][j] = mat[i][j]*mat[i][j];
			D[i][j] = vddif(coord[i],coord[j]);
                }
                Xyz[i][0] = coord[i].x;
                Xyz[i][1] = coord[i].y;
                Xyz[i][2] = coord[i].z;
        }
        for (Itno=0; Itno<50; Itno++) {
                if (!steric_xyz(D,Xyz,Dist,Strict,len,3)) break;
                new_distmat(Xyz,len,3,D);
                rmsdev_dist(D,Dist,Strict,len);
        }
        for (i=0; i<len; i++) {
                coord[i].x = Xyz[i][0];
                coord[i].y = Xyz[i][1];
                coord[i].z = Xyz[i][2];
        }
        free_matrix(Dist, len); free_matrix(Strict, len);
        free_matrix(Xyz, len);
}

float recp_area (box) Vec *box;
{
Vec	a,b,x,y;
float	ab,cd, sum, step = 0.1;
int	i, n;
	vsub(box[1],box[0],&x); vnorm(&x);
	vsub(box[3],box[2],&y); vnorm(&y);
	if (vdot(x,y) < -0.001) { printf("*NB* box edges diverge (x.y=%f)\n",vdot(x,y)); exit(1); }
	ab = vdif(box[0],box[1]);
	cd = vdif(box[2],box[3]);
	if (fabs(ab-cd) > 0.005) { printf("*NB* box has diff edges\n"); exit(1); }
	n = (int)(ab/step);
	sum = 0.0;
	vmul(&x,step); vmul(&y,step); 
	vcopy(box[0],&a); vcopy(box[2],&b);
	for (i=0; i<=n; i++)
	{ float d = vdif(a,b);
/*
		sum += 100.0/(d*d);
		sum += 10.0/d;
*/
		sum += exp(-d*d*0.01);
		vadd(a,x,&a); vadd(b,y,&b);
	}
	return sum*step;
}

float recp_end (a,b,c) Vec a,b,c;
{ /* returns inverse area for line a-b and point c */
Vec	x,y;
float	ab, sum, step = 0.1;
int	i, n;
Pv(a) NL Pv(b) NL Pv(c) NL
	vsub(b,a,&x); vnorm(&x);
	vcopy(a,&y);
	ab = vdif(a,b);
	n = (int)(ab/step);
	sum = 0.0;
	vmul(&x,step);
	for (i=0; i<=n; i++)
	{ float d = vdif(y,c);
/*
		sum += 100.0/(d*d);
		sum += 10.0/d;
*/
		sum += exp(-d*d*0.01);
		vadd(y,x,&y);
	}
	return 0.5*sum*step;
}

float	linefit(a,nat,m,w,b,den,var) Vec *a, *b, *nat; int m, w; float *den, *var;
{
int     i, j, n;
float   dot, **mom, **vec, val[3];
float	d, e, f, sum, ssum;
Vec   	c, co, cog, axis;
	mom = (float**)alloca(sizeof(float*)*3);
	vec = (float**)alloca(sizeof(float*)*3);
        for (i=0; i<3; i++) {
                mom[i] = (float*)alloca(sizeof(float)*3);
                vec[i] = (float*)alloca(sizeof(float)*3);
        }
        n = 0;
        vinit(&cog);
        for (i=m-w; i<=m+w; i++) {
		if (i < 0) continue;
                vsum(a[i],&cog);
		if (n && (vdif(a[i],a[i-1]) > BREAK)) { float d = vdif(a[i],a[i-1]);
			printf("*NB* %f BREAK between res. %d and %d\n", d,i-1,i);
			return 0.0;
		}
                n++;
        }
        vdiv(&cog,(float)n);
        for (i=0; i<3; i++) for (j=0; j<3; j++) mom[i][j] = 0.0;
        for (i=m-w; i<=m+w; i++) {
		if (i < 0) continue;
                vsub(a[i],cog,&c);
                mom[0][0] += c.x * c.x;
                mom[1][1] += c.y * c.y;
                mom[2][2] += c.z * c.z;
                mom[0][1] = mom[1][0] += c.x * c.y;
                mom[0][2] = mom[2][0] += c.x * c.z;
                mom[1][2] = mom[2][1] += c.y * c.z;
        }
	project(mom,3,val,vec);
	axis.x = vec[0][0]; axis.y = vec[0][1]; axis.z = vec[0][2];
	vadd(axis,cog, &axis);
	if (b) {
		dotOline(cog,axis,nat[m-w],b+0);
		dotOline(cog,axis,nat[m+w],b+1);
 	}
	e = sum = ssum = 0.0;
	dotOline(cog,axis,nat[m-w],&co);
	for (i=m-w; i<m+w; i++) {
		if (i < 0) continue;
		dotOline(cog,axis,nat[i+1],&c);
		d = vdif(co,c);
		f = d-e;
		sum += f;
		ssum += f*f;
		e = d;
	}
	n--;
	*den = d = sum/(float)n;
	*var = e = ssum/(float)n - d*d;
	f = val[0]/(val[1]+val[2]);
/*
	return f;
*/
	e = exp(-e*5.0);
	return e*f;
}

project (mom,Rno,val,vec)
float	**mom, **vec, *val;
int	Rno;
{
    Trimat_ Metric;
    Sqmat_ Evec;
    double *Eval, *Sqeval, *Sumeval;
    float Qu;
    double Abseval, Abssum, Newsum;
    register int i,j,Dim;
    /* build arrays */
    /* eigenvalues: real, and square roots */
    Eval=(double *) calloc(Rno,sizeof(double));
    Sqeval=(double *) calloc(Rno,sizeof(double));
    Sumeval=(double *) calloc(Rno,sizeof(double));
    /* eigenvectors: Rno*Rno */
    Metric = alloc_trimat(Rno);
    Evec=alloc_sqmat(Rno);

    /* init the actual Metric matrix */
    for (i=0; i<Rno; i++) {
        for (j=0; j<=i; j++) {
            Metric[i][j] = mom[i][j];
        }
    }

    /* get eigenvalues and eigenvectors of the metric matrix:
     in general, use the fast and precise algorithm of
     Housholder tridiag+QL transform: Numerical Recipes;
     Eigenvectors are in the rows of Evec, sorted */

    eigen_ql(Metric,Rno,Eval,Evec); /* Housholder + QL */

    /* make up Cartesian coordinates in a Dim-dimensional subspace */
    Dim = 3;
    for (j=0; j<Dim; j++) {
        if (Eval[j] > 0.0) val[j] = sqrt(Eval[j]);
		else	   val[j] = 0.0;
    }
    for (i=0; i<Rno; i++) {
        for (j=0; j<Dim; j++) {
            vec[j][i]= /* val[j]* */ Evec[j][i];     /* row eigenvectors */
	}
    }
}

pdbput (res,out,len,id,lang,temp)
Vec    *res;
FILE    *out;
int     len;
char    id;
float	*lang, *temp;
{       int     i = 0, n = 0;
        for (i=0; i<len; i++) { float a, b;
		if (lang) a = lang[i]; else a = 0.0;
		if (temp) b = temp[i]; else b = 0.0;
                fprintf(out,"ATOM%7d  CA  GLY %c%4d    %8.3f%8.3f%8.3f %5.2f %5.2f\n",
                        i+1, id, i+1, res[i].x, res[i].y, res[i].z, a, b);
        }
        fprintf(out,"TER\n");
}
