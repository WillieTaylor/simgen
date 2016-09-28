#include "util.hpp"
#include "geom.hpp"
#include "cell.hpp"
#include "data.hpp"
 
void fixSSEaxis ( Cell* cell, float weight )
{
int	level = cell->level;
int	sort = cell->sort;
int	m = cell->model;
Data	*param = Data::model+m;
Vec	x, y, now, was, axis, shift;
int     i, k, in, n = cell->kids;
float	len, fn = (float)(n-1);
float	bondCA = Data::bondCA,
	size = param->sizes[level];
int	moltype = param->moltype,
	subtype = param->subtype,
	protein = 0;
	if (n==0) return;
	if (moltype==0 && subtype==1) protein = 1;
	if (protein==1 && level==depth-1 && cell->type==2) { // protein SSE
		if (sort==2) { // beta (axis set line by i-1,i average)
			in = 0;
			axis.zero();
			was = (cell->child[0]->xyz + cell->child[1]->xyz)/2;
			for (i=2; i<n; i++) { Cell **kidi = cell->child;
				now = (kidi[i-1]->xyz + kidi[i]->xyz)/2;
				axis += now - was;
				was = now;
				in++;
			}
			if (in<1) return;
		}
		if (sort==1) { // alpha (axis line set by i-1,i,i+1 average)
			in = 0;
			axis.zero();
			was = (cell->child[0]->xyz + cell->child[1]->xyz + cell->child[2]->xyz)/3;
			for (i=2; i<n-1; i++) { Cell **kidi = cell->child;
				now = (kidi[i-1]->xyz + kidi[i]->xyz + kidi[i+1]->xyz)/3;
				axis += now - was;
				was = now;
				in++;
			}
			if (in<1) return;
		}
		if (sort==0) { // loop (set in direction of first half to last half)
			if (n<2) return;
			m = n/2;
			x.zero(); k = 0;
			for (i=0; i<m; i++) { x += cell->child[i]->xyz; k++; }
			x /= (float)k;
			y.zero(); k = 0;
			for (i=m; i<n; i++) { y += cell->child[i]->xyz; k++; }
			y /= (float)k;
			axis = y-x;
			weight *= 0.1;
		}
		// reset the end-points relative to the centre
		if (sort==0) len = loopAXIS*size*sqrt(fn);
		if (sort==1) len = 0.7*alphAXIS*fn*bondCA/3.8;
		if (sort==2) len = betaAXIS*fn*bondCA/3.8;
		axis.setVec(len*0.5);
		x = cell->xyz - axis;
		shift = x - cell->endN;
		cell->endN += shift*weight;
		y = cell->xyz + axis;
		shift = y - cell->endC;
		cell->endC += shift*weight;
		cell->ends = 1;
		if (sort) { // refine alpha and beta children using cent values
			FOR(i,n) { Cell *kidi = cell->child[i]; Vec dif; float rms1, rms2;
				dif.x = kidi->xyz|cell->endN;
				dif.y = kidi->xyz|cell->xyz;
				dif.z = kidi->xyz|cell->endC;
				dif -= kidi->cent; // difference between have and want (cent) 
				rms1 = sqrt(dif*dif);
				shift.zero();
				shift += (cell->endN - kidi->xyz)*dif.x;
				shift += (cell->xyz  - kidi->xyz)*dif.y;
				shift += (cell->endC - kidi->xyz)*dif.x;
				shift *= weight/3;
				now = kidi->xyz + shift;
				dif.x = now|cell->endN;
				dif.y = now|cell->xyz;
				dif.z = now|cell->endC;
				dif -= kidi->cent;
				rms2 = sqrt(dif*dif);
				if (rms1-rms2 > NOISE) { // apply the shift
					kidi->xyz = now;
				}
			}
		} 
	}
	FOR(i,n) fixSSEaxis(cell->child[i], weight);
}
 
void fixRNAaxis ( Cell* cell, float weight )
{
int	level = cell->level;
int	type = cell->type;
int	sort = cell->sort;
int	m = cell->model;
Data	*param = Data::model+m;
int     n = cell->kids;
int	moltype = param->moltype,
	subtype = param->subtype;
	if (n==0) return;
	if (moltype==1 && sort==1 && level==depth-2 && type==2) // nucleic SSE
        { int   to,in = -1; float fn = 0.0;
	  Vec	shift, axis, cent, endN, endC;
		// stem axis set by end pairs of basepairs
                cent.zero();
                FOR(i,n) { Cell *b = cell->child[i];
                	if (b->sort == 0) continue; // use only basepairs
                        cent += b->xyz; fn += 1.0;
                        if (in<0) in = i;
                        to = i;
                }
                cell->xyz = cent/fn;
                endN = cell->child[in]->xyz & cell->child[in+1]->xyz;
                endC = cell->child[to]->xyz & cell->child[to-1]->xyz;
                axis = endC - endN;            // set axis direction
		fn = 0.23 * (fn-1.0);	// RNA = 2.3 rise/bp (1/10 for scale)
		if (fn < 1.0) fn = 1.0;
		axis.setVec(fn*0.5);	// set to half length
		// shift endN towards target
		shift = (cell->xyz - axis) - cell->endN;
		cell->endN += shift*weight;
		// shift endC towards target
		shift = (cell->xyz + axis) - cell->endC;
		cell->endC += shift*weight;
	}
	FOR(i,n) fixRNAaxis(cell->child[i], weight);
}
 
 
void fixBPbonds ( Cell***, int, int, float );
void fixBPbumps ( Cell***, int, int, float );

void fixRNAstem ( Cell* cell, float weight )
{
float	thgiew = 1.0 - weight;
int	level = cell->level;
int	sort = cell->sort;
int	mod = cell->model;
Data	*param = Data::model+mod;
int     m,n, kids = cell->kids;
float	size = param->sizes[level];
int	moltype = param->moltype,
	subtype = param->subtype;
Vec	w,c, mid, axis, shift, spoke;
float	BPrise, PPrise, radius, twist;
float	d, d0,dN;
Seg	pole;
int	fix = 1;
	if (kids==0) return;
	radius = 1.5;
	BPrise = 0.55;
	PPrise = 0.51;// 0.40;
	if ( drand48() < 0.5 // don't refine everyone at the same time
		&& moltype==1 && cell->type==2 && cell->sort && level==depth-2) // RNA stem
	{ Cell  ***stem = new Cell**[kids];	// make a temp copy in stem
		FOR(i,kids) stem[i] = new Cell*[2];
		n = 0;
		FOR(i,kids) { Cell *pi = cell->child[i];
			// discard loop segments
			if (pi->sort == 0) continue;
			// copy cell coordinates to stem
			stem[n][0] = new Cell;
			stem[n][0]->xyz = pi->child[0]->xyz;
			stem[n][1] = new Cell;
			stem[n][1]->xyz = pi->child[1]->xyz;
			n++;
		}
		fixBPbumps(stem,n,fix,weight);
		fixBPbonds(stem,n,fix,weight);
		m = n-1;
		d0 = (stem[0][0]->xyz & stem[0][1]->xyz)  | cell->endN;
		dN = (stem[m][0]->xyz & stem[m][1]->xyz)  | cell->endN;
		if (d0 < dN) pole = Seg(cell->endN,cell->endC);
		if (d0 > dN) pole = Seg(cell->endC,cell->endN);
		axis = (pole.B-pole.A).getNorm();		// unit vector along the axis
		// move the Ps to the tube surface and reset angle between paired bases
		twist = 1.25; //1.0; //1.25;		// angle between basepaired Ps
		FOR(i,n)
		{ Cell *ci = stem[i][0], *wi = stem[i][1];
		  Vec	sc, sw;
int k=0;
			// get the image of the Ps on the axis for each strand
			c = ci->xyz.vec_on_line(pole);	// image of c on the axis
			w = wi->xyz.vec_on_line(pole);	// image of w on the axis
{ Vec v =      c; printf("ATOM%7d  CA  GLY C%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,k++,v.x,v.y,v.z); }
{ Vec v =      w; printf("ATOM%7d  CA  GLY C%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,k++,v.x,v.y,v.z); }
			mid = c & w;
			sc = ((ci->xyz - c)^axis).getVec(radius);	// 90 to c spoke
			sw = ((w - wi->xyz)^axis).getVec(radius);	// 90 to w spoke (on same side)
Pr(sc*sw) NL
			spoke = mid + (sc+sw).getVec(radius);	 // bisector of WC angle
			c = spoke.get_rot(pole,-twist) + axis*BPrise*0.5; 
			stem[i][0]->xyz = ci->xyz & c;
			w = spoke.get_rot(pole, twist) - axis*BPrise*0.5;
			stem[i][1]->xyz = wi->xyz & w;
{ Vec v =      c; printf("ATOM%7d  CA  GLY C%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,k++,v.x,v.y,v.z); }
{ Vec v =      w; printf("ATOM%7d  CA  GLY C%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,k++,v.x,v.y,v.z); }
{ Vec v =    mid; printf("ATOM%7d  CA  GLY C%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,k++,v.x,v.y,v.z); }
{ Vec v = sc+mid; printf("ATOM%7d  CA  GLY C%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,k++,v.x,v.y,v.z); }
{ Vec v = sc+mid; printf("ATOM%7d  CA  GLY C%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,k++,v.x,v.y,v.z); }
		}
/*
		fixBPbumps(stem,n,fix,weight);
		fixBPbonds(stem,n,fix,weight);
		// set the twist
		twist = -0.6;//-0.55;
		FOR(i,n) { Cell *ci,*wi, *cj,*wj;
			if (drand48()<0.5) { // fix up
				if (i > n-2) continue;
				ci = stem[i][0];   wi = stem[i][1];
				cj = stem[i+1][0]; wj = stem[i+1][1];
			} else {	// fix down
				if (i < 1) continue;
				ci = stem[i-1][0]; wi = stem[i-1][1];
				cj = stem[i][0];   wj = stem[i][1];
			}
			c = ci->xyz.get_rot(pole,twist);
			c += axis*PPrise;
			cj->xyz = cj->xyz & c;
			w = wi->xyz.get_rot(pole,twist);
			w += axis*PPrise;
			wj->xyz = wj->xyz & w;
		}
		fixBPbumps(stem,n,fix,weight);
		fixBPbonds(stem,n,fix,weight);
*/
printf("TER  GLY\n");
sleep(1);
{ Vec wi = cell->endN; int i = 10;
printf("ATOM%7d  CA  GLY C%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,wi.x,wi.y,wi.z);
}
{ Vec wi = cell->endC; int i = 20;
printf("ATOM%7d  CA  GLY C%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,wi.x,wi.y,wi.z);
}
printf("TER  GLY\n");
FOR(i,cell->kids) { Cell *pi = cell->child[i], *ci = pi->child[0];
printf("ATOM%7d  CA  GLY A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,ci->xyz.x,ci->xyz.y,ci->xyz.z);
}
FOR(i,cell->kids) { Cell *pi = cell->child[i], *ci = pi->child[1];
printf("ATOM%7d  CA  GLY A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,ci->xyz.x,ci->xyz.y,ci->xyz.z);
}
printf("TER  GLY\n");
FOR(i,n) { Cell *ci = stem[i][0];
printf("ATOM%7d  CA  GLY B%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i+1,ci->xyz.x,ci->xyz.y,ci->xyz.z);
}
FOR(i,n) { Cell *ci = stem[i][1];
printf("ATOM%7d  CA  GLY B%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i+1,ci->xyz.x,ci->xyz.y,ci->xyz.z);
}
printf("TER  GLY\n");
exit(1);
		// copy back (weighted) basepair positions
		n = 0;
		FOR(i,kids) { Cell *pi = cell->child[i];
			if (pi->sort == 0) continue;
			pi->child[0]->xyz = stem[n][0]->xyz*weight + pi->child[0]->xyz*thgiew;
			pi->child[1]->xyz = stem[n][1]->xyz*weight + pi->child[1]->xyz*thgiew;
			n++;
		}
	}
	FOR(i,kids) fixRNAstem(cell->child[i], weight);
}

void fixBPbonds ( Cell ***stem, int n, int fixes, float kick ) {
// patch-up local BP and backbone geometry
float	bondBP  = 2.5;
float	bondPP  = 0.9;
float	bondPxP = 1.7;
//float	kick = 0.1;
	FOR(j,fixes) { // fixing cycles
		FOR(i,n) // fix BP separation
		{ Cell *ci = stem[i][0], *wi = stem[i][1];
			part2cells(ci,wi,bondBP,kick);
		}
		DO1(i,n-2) // fix PxP separation on Watson strand
		{ Cell *bi = stem[i-1][0], *ci = stem[i+1][0];
			part2cells(bi,ci,bondPxP,kick);
		}
		DO1(i,n-1) // fix PP separation on Watson strand
		{ Cell *bi = stem[i-1][0], *ci = stem[i][0];
			part2cells(bi,ci,bondPP,kick);
		}
		DO1(i,n-2) // fix PxP separation on Crick strand
		{ Cell *bi = stem[i-1][1], *ci = stem[i+1][1];
			part2cells(bi,ci,bondPxP,kick);
		}
		DO1(i,n-1) // fix PP separation on Crick strand
		{ Cell *bi = stem[i-1][1], *ci = stem[i][1];
			part2cells(bi,ci,bondPP,kick);
		}
	}
}

void fixBPbumps ( Cell ***stem, int n, int fixes, float kick ) {
// patch-up local BP and backbone geometry
float	d, bump = 2.0;
//float	kick = 0.1;
	FOR(k,fixes) { // fixing cycles
		FOR(i,n) // fix BP separation
		{ Cell *ci = stem[i][0], *wi = stem[i][1];
			FOR(j,n) // fix bumps with Watson i
			{ Cell  *cj = stem[j][0], *wj = stem[j][1];
				d = ci->xyz|cj->xyz;
				if (d<bump && (i-j>3 || j-i>3)) { // C-C strand bump
					part2cells(ci,cj,bump,kick);
				}
				d = ci->xyz|wj->xyz;
				if (d<bump) {			// C-W strand bump
					part2cells(ci,wj,bump,kick);
				}
				d = wi->xyz|wj->xyz;
				if (d<bump && (i-j>3 || j-i>3)) { // W-W strand bump
					part2cells(wi,wj,bump,kick);
				}
				d = wi->xyz|cj->xyz;
				if (d<bump) {			// W-C strand bump
					part2cells(wi,cj,bump,kick);
				}
			}
		}
	}
}
/*
		new1 = (Vec*)alloca(sizeof(Vec)*cell->kids);
		new2 = (Vec*)alloca(sizeof(Vec)*cell->kids);
//for (i=0; i<cell->kids; i++) { Cells *ci = cell->child[i];
//printf("ATOM%7d  CA  GLY A%4d     %7.3f %7.3f %7.3f  0.00 %5.2f\n", i,i,ci->xyz.x,ci->xyz.y,ci->xyz.z,(float)ci->lost);
//}
//printf("TER ATOM\n");
		vsub(cell->endC,cell->endN,&axis);
//{ float d = vdif(cell->child[35]->xyz,cell->child[18]->xyz); Pr(d) NL }
//return;
		for (i=0; i<n; i++) { Cells *cat = cell->child[i], *cdt, *cj;
			if (cat->link==0) continue;
			cdt = cat->link[0];
			if (cdt==0) continue;
			for (j=cdt->id-5; j<=cdt->id; j++) { float d;
				if (j<0) continue;
				if (abs(cat->id-j) < 3) continue;
				cj = cell->child[j];
				d = vdif(cat->xyz,cj->xyz);
				if (d<clash) part2cells(cat,cj,clash,0.1,0);
//Pi(cat->id) Pi(cdt->id) Pi(cj->id) Pr(d) NL
			}
			for (j=cat->id+5; j>=cat->id; j--) { float d;
				if (j>=n) continue;
				if (abs(cdt->id-j) < 3) continue;
				cj = cell->child[j];
				d = vdif(cdt->xyz,cj->xyz);
				if (d<clash) part2cells(cdt,cj,clash,0.1,0);
//Pi(cdt->id) Pi(cat->id) Pi(cj->id) Pr(d) NL
			}
		}
//return;
		for (i=0; i<n; i++)
		{ Cells *cat = cell->child[i], *cdt;
		  Vec	a,b,c,d;
			if (cat->link==0) continue;
			cdt = cat->link[0];
			if (cdt==0) continue;
			vcopy(cat->xyz,&a); vcopy(cdt->xyz,&d);
			ab = dotOline(cell->endC,cell->endN,a,&b);
			cd = dotOline(cell->endC,cell->endN,d,&c);
/////////
if (i) {
float dis = vdif(a,cell->child[i-1]->xyz),
bond = 0.1*(float)(sizes[level+1]+bonds[level+1]);
bc = vdif(b,c);
ad = vdif(a,d);
Pr(dis) Pr(bond) Pr(radius) Pr(ab) Pr(cd) Pr(bc) Pr(ad) NL
}
////////
			vave(b,c,&p); vave(a,d,&q); vsub(q,p,&r);	// r = radial spoke
			vnorm(&r); vmul(&r,radius);
			rise = 1.00/2.0;
			vsub(c,b,&w); vnorm(&w); vmul(&w,rise);	// w = half axial rise
			vsub(p,w,&v); vadd(r,v,&at);	// at is on same level as cat
			vadd(p,w,&u); vadd(r,u,&dt);	// dt is on same level as cdt
			twist = 1.70/2.0;
			rotate(b,c,&at,-twist);
			rotate(b,c,&dt, twist);
			vave(a,at,new1+i); // new1 = average old and ideal at i
			vave(d,dt,new2+i); // new2 = average at linked from i
		}
		for (i=0; i<n; i++) { Cells *cat = cell->child[i], *cdt;
			if (cat->link==0) continue;
			cdt = cat->link[0];
			if (cdt==0) continue;
			vsub(new1[i],cat->xyz,&shift); vmul(&shift,weight);
			vadd(cat->xyz,shift,&(cat->xyz));
			vsub(new2[i],cdt->xyz,&shift); vmul(&shift,weight);
			vadd(cdt->xyz,shift,&(cdt->xyz));
		}
//for (i=0; i<cell->kids; i++) { Cells *ci = cell->child[i];
//printf("ATOM%7d  CA  ALA A%4d     %7.3f %7.3f %7.3f  0.00 %5.2f\n", i,i,ci->xyz.x,ci->xyz.y,ci->xyz.z,(float)ci->lost);
//}
//exit(1);
	}
	for (i=0; i<n; i++) fixRNAstem(cell->child[i], level+1, weight);
}
 
void fixRNAaxis ( Cells* cell, int level, float weight )
{
Vec	x, mid, axis, axwas, shift, newN, newC;
int     i,j,k, m, n = cell->kids;
float	len;
	m = cell->model*M*N+N; sizes = data+m+N*3;
	moltype = data[m];
	if (n==0) return;
	if (moltype==1 && cell->type==2) {
		if (cell->sort==0) { // loop: set axis by N-C half location
			m = n/2;
			k = 0;
			vinit(&mid);
			for (i=0; i<m; i++) {
				vsum(cell->child[i]->xyz, &mid); k++;
			}
			vdiv(&mid,(float)k);
			vcopy(mid,&newN);
			k = 0;
			vinit(&mid);
			for (i=m; i<n; i++) {
				vsum(cell->child[i]->xyz, &mid); k++;
			}
			vdiv(&mid,(float)k);
			vcopy(mid,&newC);
		}
		if (cell->sort==1) { Vec *mids; int *use;
			mids = (Vec*)alloca(sizeof(Vec)*n);
			use  = (int*)alloca(sizeof(int)*n);
			FOR(i,n) {
				cell->child[i]->lost = 1;
				mids[i].x = 999.9;
			}
			FOR(i,n) { Cells *ci = cell->child[i];
				if (ci->link==0 || ci->link[0]==0 ) continue; // no links
				if (ci->link[0]->parent != cell) continue; // linked outside
				ci->lost = 0; ci->link[0]->lost = 0;
			}
			for (j=0; j<n; j++) { Cells *ci, *cj, *ck; int skip;
				i = j-6; k = j+6; // +/-6 = half turn up/down
				if (i<0 || k>=n) continue;
				ci = cell->child[i]; cj = cell->child[j]; ck = cell->child[k];
				if (ci->lost || cj->lost || ck->lost) continue;
				skip = 0;
				for (m=i; m<=j; m++) { if (cell->child[m]->lost) skip = 1; }
				for (m=j; m<=k; m++) { if (cell->child[m]->lost) skip = 1; }
				if (skip) continue;
				vave(ci->xyz,ck->xyz,&mid);
				vave(mid,cj->xyz,mids+j);
//printf("ATOM%7d  CA  VAL A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,mids[j].x,mids[j].y,mids[j].z);
			}
//printf("TER ATOM\n");
			vcopy(cell->endN,&newN);
			vcopy(cell->endC,&newC);
			k = min(3,n/4);
			FOR(i,k) { float d, dmin;
				// shift newN towards nearest remaining mid-point
				dmin = 999.9;
				FOR(j,n) {
					if (mids[j].x>999.0) continue;	
					d = vdif(mids[j],newN);
					if (d<dmin) { dmin = d; m = j; }
				}
				vave(mids[m],newN,&newN);
//printf("ATOM%7d  CA  LYS A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i+1,newN.x,newN.y,newN.z);
//printf("TER ATOM\n");
				mids[m].x = 999.9;
				// shift newC towards nearest remaining mid-point
				dmin = 999.9;
				FOR(j,n) {
					if (mids[j].x>999.0) continue;	
					d = vdif(mids[j],newC);
					if (d<dmin) { dmin = d; m = j; }
				}
				vave(mids[m],newC,&newC);
//printf("ATOM%7d  CA  ASP A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,-i-1,newC.x,newC.y,newC.z);
//printf("TER ATOM\n");
				mids[m].x = 999.9;
			}
		}
//printf("ATOM%7d  CA  ALA A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,cell->endN.x,cell->endN.y,cell->endN.z);
//printf("ATOM%7d  CA  ALA A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,cell->endC.x,cell->endC.y,cell->endC.z);
//printf("TER ATOM\n");
//weight=1.0;
		vsub(newC,newN,&axis);
		if (cell->sort==0) len = 0.2*sqrt((float)n);
		// RNA = 2.3 rise/bp --> 0.115 (1/2 for duplex, 1/10 for scale)
		if (cell->sort==1) len = 0.115*(float)n;
		vnorm(&axis); vmul(&axis,len*0.5);
		vadd(cell->xyz,axis,&x); // x is near C
//printf("ATOM%7d  CA  SER A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,x.x,x.y,x.z);
//printf("TER ATOM\n");
		vsub(x,cell->endC, &shift); vmul(&shift,weight);
		vsum(shift, &(cell->endC));
		vsub(cell->xyz,axis,&x); // x is near N
//printf("ATOM%7d  CA  HIS A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,x.x,x.y,x.z);
//printf("TER ATOM\n");
		vsub(x,cell->endN, &shift); vmul(&shift,weight);
		vsum(shift, &(cell->endN));
		cell->ends = 1;
		if (cell->sort==1) { int bad, bent; // fix rise/bp
//for (i=0; i<cell->kids; i++) { Cells *ci = cell->child[i];
//printf("ATOM%7d  CA  GLY A%4d     %7.3f %7.3f %7.3f  0.00 %5.2f\n", i,i,ci->xyz.x,ci->xyz.y,ci->xyz.z,(float)ci->lost);
//} i=0;
//printf("TER ATOM\n");
//printf("ATOM%7d  CA  CYS A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,cell->endN.x,cell->endN.y,cell->endN.z);
//printf("ATOM%7d  CA  CYS A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,cell->endC.x,cell->endC.y,cell->endC.z);
//exit(1);
			bad = bent = 0;
			m = n/2; // assume half is up, half is down 
			vnorm(&axis);
			for (i=1; i<n-1; i++)
			{ Cells *a, *b, *c;
			  Vec	ax,bx,cx, da,db,dc, mid;
	  		  float ab,ac,bc, d, an,cn, rise = 0.30;
			  int	swap = 0;
				a = cell->child[i-1];
				b = cell->child[i  ];
				c = cell->child[i+1];
				if (a->lost || b->lost || c->lost) continue; // skip loop (no bp)
				d = vdif(a->xyz,c->xyz);
				if (i>1 && i<n-2 && d<1.0) bent++;
				dotOline(cell->endC,cell->endN,a->xyz,&ax);
				dotOline(cell->endC,cell->endN,b->xyz,&bx);
				dotOline(cell->endC,cell->endN,c->xyz,&cx);
				ab=vdif(ax,bx); bc=vdif(bx,cx); ac=vdif(ax,cx); d=ac-ab-bc;
				if (d > -NOISE && ab>0.1 && bc>0.1) continue;
				// bases are out of order on the axis or too close (rise/3)
				an = vdif(ax,cell->endN);
				cn = vdif(cx,cell->endN);
				d = an-cn;
				if (b->id < m) { // ascending (d = an-cn = +ve)
					if (d<0.0) swap = 1;
				} else {	 // decending (d = an-cn = -ve)
					if (d>0.0) swap = 1;
				}
				if (swap) bad++;
				vave(ax,cx, &mid);
				vsub(ax,mid, &da); vnorm(&da); vmul(&da,rise);
				vsub(cx,mid, &dc); vnorm(&dc); vmul(&dc,rise);
				if (swap) { vcopy(da,&x); vcopy(dc,&da); vcopy(x,&dc); }
				vsum(mid,&da);
				vsum(mid,&dc);
				vsub(da, ax, &da); // da = shift to new a
				vsub(mid,bx, &db); // db = shift to centre b
				vsub(dc, cx, &dc); // dc = shift to new c
				vmul(&da,weight); vmul(&db,weight); vmul(&dc,weight);
				vsum(da,&(a->xyz));
				vsum(db,&(b->xyz));
				vsum(dc,&(c->xyz));
			}
if (bad+bent>10) { Pi(bad) Pi(bent) Pi(cell->parent->id) NL }
			FOR(j,bad+bent) repair(cell);
		}
		FOR(i,n) cell->child[i]->lost = 0;
	}
	for (i=0; i<n; i++) fixRNAaxis(cell->child[i], level+1, weight);
}

repair ( Cells *cell )
{
float	a,b,c, AB=2.8, C=3.7;
Vec	x,y,z, r,t;
float	rad = 1.0, bond = 0.56, turn = 3.76, twist = 0.0;
int	i,j,k,m,n, jj, len = cell->kids, mid = len/2;
int	rank[99];
float	skew[99];
Cells *ci,*cj,*ck,*cm,*cn;
Vec	vi,vj,vk,vm,vn;
Vec	ui,uj,uk,um,un;
Vec	xj, yj;
//////////////////
Pi(mid) NL
i = 0;
printf("ATOM%7d  CA  ARG A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,cell->endN.x,cell->endN.y,cell->endN.z);
i = 1;
printf("ATOM%7d  CA  GLU A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,cell->endC.x,cell->endC.y,cell->endC.z);
printf("TER ATOM\n");
FOR(i,len) { Cells *ci = cell->child[i];
printf("ATOM%7d  CA  CYS A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,ci->id,ci->xyz.x,ci->xyz.y,ci->xyz.z);
}
printf("TER ATOM\n");
//////////////////
//rad = 0.76;
	vsub(cell->endC,cell->endN,&x); vnorm(&x);
	vcopy(x,&r); vmul(&r,turn/11.0); // rise/base
	vcopy(x,&t); vmul(&t,turn*0.5); // half turn
	FOR(j,len) {
		skew[j] = -999.9;
		i = j-1;  k = j+1;
		m = j-6; n = j+6;
		if (m < 0) continue;
		if (n >= len) continue;
		if (m<mid && n>mid) continue;
		ci = cell->child[i]; vcopy(ci->xyz,&vi);
		cj = cell->child[j]; vcopy(cj->xyz,&vj);
		ck = cell->child[j]; vcopy(ck->xyz,&vk);
		cm = cell->child[m]; vcopy(cm->xyz,&vm);
		cn = cell->child[n]; vcopy(cn->xyz,&vn);
//a = vdif(vj,vm); b = vdif(vj,vn); c = vdif(vn,vm);
//Pi(j) Pr(a) Pr(b) Pr(c) NL
		a = vdif(vj,vm); a = a-AB; a*=a;
		b = vdif(vj,vn); b = b-AB; b*=b;
		c = vdif(vn,vm); c = c-C;  c*=c;
		skew[j] = -a-b-c; // -ve sos of deviations from ideal triangle
	}
	sort(0,skew,0,rank,len,1);
	FOR(jj,n) { int new;
		j = rank[jj];
		if (skew[j] < -1.0) break;
		i = j-1;  k = j+1;
		m = j-6; n = j+6;
		if (m < 0) continue;
		if (n >= len) continue;
		if (m<mid && n>mid) continue;
		ci = cell->child[i]; vcopy(ci->xyz,&vi);
		cj = cell->child[j]; vcopy(cj->xyz,&vj);
		ck = cell->child[k]; vcopy(ck->xyz,&vk);
		cm = cell->child[m]; vcopy(cm->xyz,&vm);
		cn = cell->child[n]; vcopy(cn->xyz,&vn);
		dotOline(cell->endC,cell->endN,vj,&xj);
		vsub(vj,xj,&y); vnorm(&y); vmul(&y,rad);
		vprod(x,y, &z); vnorm(&z); vmul(&z,bond); 
		vadd(xj,y,&uj); // new j
		vsub(xj,y,&yj); // point opposite j
if (1) {
//printf("ATOM%7d  CA  ALA A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", j,j,uj.x,uj.y,uj.z);
		if (j>mid) {
			vsub(uj,z,&ui); vsub(ui,r,&ui);
			vadd(uj,z,&uk); vadd(uk,r,&uk);
		} else {
			vsub(uj,z,&uk); vsub(uk,r,&uk);
			vadd(uj,z,&ui); vadd(ui,r,&ui);
		}
		vave(vi,ui,&(ci->xyz));
		vave(vk,uk,&(ck->xyz));
//printf("ATOM%7d  CB  HIS A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", j,-i,ui.x,ui.y,ui.z);
//printf("ATOM%7d  CB  ASP A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", j,-k,uk.x,uk.y,uk.z);
		twist = 0.01;
		new = 0;
		if (j<mid && m<mid) { vadd(yj,t,&um); new = 1; }
		if (j>mid && m>mid) { vsub(yj,t,&um); new = 1; }
		if (new) {
			rotate(cell->endC,cell->endN,&um,twist);
			vave(vm,um,&(cm->xyz));
		}
		new = 0;
		if (j>mid && n>mid) { vadd(yj,t,&un); new = 1; }
		if (j<mid && n<mid) { vsub(yj,t,&un); new = 1; }
		if (new) {
			rotate(cell->endC,cell->endN,&un,twist);
			vave(vn,un,&(cn->xyz));
		}
//printf("ATOM%7d  CB  HIS A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", j,-m,um.x,um.y,um.z);
//printf("ATOM%7d  CB  ASP A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", j,-n,un.x,un.y,un.z);
//exit(1);
} else {
//printf("ATOM%7d  CA  VAL A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", j,j,uj.x,uj.y,uj.z);
}
	}
}

oldrepair ( Cells *cell )
{
Vec	x,y,z, r,t;
float	rad = 0.78, bond = 0.59, turn = 2.3, rise = 0.2;
int	i,j,k,m,n, len = cell->kids, mid = len/2;
FOR(i,len) { Cells *ci = cell->child[i];
printf("ATOM%7d  CA  CYS A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,ci->xyz.x,ci->xyz.y,ci->xyz.z);
}
printf("TER ATOM\n");
	vsub(cell->endC,cell->endN,&x); vnorm(&x);
	vcopy(x,&r); vmul(&r,rise);
	vcopy(x,&t); vmul(&t,turn);
	FOR(j,len)
	{ Cells *ci,*cj,*ck,*cm,*cn;
	  Vec	vi,vj,vk,vm,vn;
	  Vec	ui,uj,uk,um,un;
	  Vec	xj;
		i = j-1;  k = j+1;
		m = j-10; n = j+10;
		if (m < 0) continue;
		if (n >= len) continue;
		ci = cell->child[i]; vcopy(ci->xyz,&vi);
		cj = cell->child[j]; vcopy(cj->xyz,&vj);
		ck = cell->child[j]; vcopy(ck->xyz,&vk);
		cm = cell->child[m]; vcopy(cm->xyz,&vm);
		cn = cell->child[n]; vcopy(cn->xyz,&vn);
		dotOline(cell->endC,cell->endN,vj,&xj);
		vsub(vj,xj,&y); vnorm(&y); vmul(&y,rad);
		vprod(x,y, &z); vnorm(&z); vmul(&z,bond); 
		vadd(xj,y,&uj); 
printf("ATOM%7d  CA  VAL A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,uj.x,uj.y,uj.z);
//		vsub(uj,z,&ui);
//		vadd(uj,z,&uk);
///////////////////////
Pr(vdif(vi,ui));
Pr(vdif(vj,uj));
Pr(vdif(vk,uk));
NL
		if (j<mid) {
			vadd(ui,r,&ui);
			vsub(uk,r,&uk);
		} else {
			vsub(ui,r,&ui);
			vadd(uk,r,&uk);
		}
///////////////////////
//		vcopy(vm, &um); vcopy(vn, &un);
//		if (j<mid && m<mid) vadd(uj,t,&um);
//		if (j>mid && n>mid) vadd(uj,t,&un);
		//if (j<mid && m<mid) vsub(uj,t,&um);
		//if (j>mid && n>mid) vsub(uj,t,&un);
Pi(j)
Pr(vdif(vi,ui));
Pr(vdif(vj,uj));
Pr(vdif(vk,uk));
Pr(vdif(vm,um));
Pr(vdif(vn,un));
NL
//		vave(ui,vi,&(ci->xyz));
		vave(uj,vj,&(cj->xyz));
//		vave(uk,vk,&(ck->xyz));
//		vave(um,vm,&(cm->xyz));
//		vave(un,vn,&(cn->xyz));
	}
printf("TER ATOM\n");
FOR(i,len) { Cells *ci = cell->child[i];
printf("ATOM%7d  CA  ALA A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,ci->xyz.x,ci->xyz.y,ci->xyz.z);
}
exit(1);
NL
}


remake ( Cells *cell )
{
int	i, n = cell->kids, m = n/2;
Vec	axis, step, new;
float	x,y,z, c,s,t, sign = 1.0, d = 1.0;
int	along;
Pt(remake) NL
	vcopy(cell->endN, &axis);
	vsub(cell->endC, cell->endN, &step); vdiv(&step,(float)m);
	x = step.x; x*=x; y = step.y; y*=y; z = step.z; z*=z;
	if (x>y && x>z) { along = 1; if (step.x<0.0) sign = -1.0; }
	if (y>x && y>z) { along = 2; if (step.y<0.0) sign = -1.0; }
	if (z>x && z>y) { along = 3; if (step.z<0.0) sign = -1.0; }
//Pt(ATOM) Pi(along) Pr(sign) NL
	FOR(i,n) { Cells *ci = cell->child[i];
		t = 0.5*(float)i;
		s = d*sin(t); c = d*cos(t); 
		if (i<m) {
			vadd(axis,step, &axis); 
		} else {
			vsub(axis,step, &axis);
			t=c; c=s; s=t;
		}
		s *= sign;
		vcopy(axis,&new);
		if (along==1) { new.y += c; new.z += s; }
		if (along==2) { new.z += c; new.x += s; }
		if (along==3) { new.x += c; new.y += s; }
//printf("ATOM%7d  CA  VAL A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,ci->xyz.x,ci->xyz.y,ci->xyz.z);
		vave(new,ci->xyz, &(ci->xyz));
	}
//exit(1);
}
 
void fixDOMaxis ( Cells* cell, int level, float weight )
{
Vec	x, old, new, shift;
int     i, m, n = cell->kids;
float	len;
	m = cell->model*M*N+N; sizes = data+m+N*3;
	len = 0.1*(float)sizes[level];
	if (n==0) return;
	if (cell->type==3 && cell->child[0]->type==2) { // reset ellipsoid axis if there are SSEs inside
		vinit(&new);
		vsub(cell->endC,cell->endN, &old);	// current N--C axis
		if (cell->bond) {	// pick based on neighbours if none given
			if (vmod(old) < NOISE) vsub(cell->bro->xyz,cell->sis->xyz, &old);
			if (vmod(old) < NOISE) vsub(cell->child[0]->xyz,cell->xyz, &old);
		} else {		// otherwise random
			vrset(&old,1.0);
		}
		// recalculate axis as average SSE axis
		for (i=0; i<n; i++) { Cells *child = cell->child[i];
			if (child->sort == 0) continue;	// exclude loops
			vsub(child->endC,child->endN, &x);
			if (vdot(x,old)>0.0) {
				vadd(new, x, &new);	// sum parallel SSE
			} else {
				vsub(new, x, &new);	// sum antipar. SSE
			}
		}
		len *= 0.1*(float)(cell->sort+10-E/2); // for E ellipsoid sorts
		vnorm(&new); vmul(&new,len*0.5);
		vsub(cell->xyz, new, &x);
		vsub(x,cell->endN, &shift); vmul(&shift,weight);
		vsum(shift, &(cell->endN));
		vadd(cell->xyz, new, &x);
		vsub(x,cell->endC, &shift); vmul(&shift,weight);
		vsum(shift, &(cell->endC));
		separate(&(cell->endN),&(cell->endC),len,weight);
		cell->ends = 1;
///////////////////////
	} else {
		if (vdif(cell->endN,cell->endC) > NOISE) { // endpoints exist
			cell->ends = 1;
		}	
///////////////////////
	}
	for (i=0; i<n; i++) fixDOMaxis(cell->child[i], level+1, weight);
}
 
void reset2axis ( Cells* cell, int level, float weight )
{
int     i, m, n, type, sort;
Vec	new, newN, newC, shift;
float	size, len;
	n = cell->kids;
	if (n==0) return;
	if (n==1) return;
	m = cell->model*M*N+N; sizes = data+m+N*3;
	moltype = data[m];
	type = abs(cell->type);
	sort = cell->sort;
	if (level<depth-1 && cell->child[0]->cent.x > 0.0) { // reset positions to axes (using cent distances)
		FOR(i,n)
		{ Cells *kidi = cell->child[i];
		  float a = kidi->cent.x, b = kidi->cent.y, c = kidi->cent.z;
		  Vec	A, B, C;
			if (b<NOISE) continue;
			vatA(cell->endN,kidi->xyz,&A,a);
			vatA(cell->xyz ,kidi->xyz,&B,b);
			vatA(cell->endC,kidi->xyz,&C,c);
			vinit(&new);
			vsum(A,&new); vsum(B,&new); vsum(C,&new);
			vdiv(&new,3.0);
			vsub(new,kidi->xyz, &shift); vmul(&shift,weight);
			moveCell(kidi,shift,1);
		}
	}
	FOR(i,n) reset2axis(cell->child[i], level+1, weight);
}

int goodTri ( float a, float b, float c )
{
	if (b+c < a) return 0;
	if (c+a < b) return 0;
	if (a+b < c) return 0;
	if (a>b) { if (a-b > c) return 0; } else { if (b-a > c) return 0; }
	if (b>c) { if (b-c > a) return 0; } else { if (c-b > a) return 0; }
	if (c>a) { if (c-a > b) return 0; } else { if (a-c > b) return 0; }
	return 1;
}
 
void fixANYaxis ( Cells* cell, int level, float weight )
{
int     i, m, n, inN, inC, set, type, sort;
Vec	new, newN, newC, shift;
float	size, len;
	n = cell->kids;
	if (n==0) return;
	m = cell->model*M*N+N; sizes = data+m+N*3;
	moltype = data[m];
	type = abs(cell->type);
	sort = cell->sort;
	size = 0.1*(float)sizes[level];
	if (type==1) { // sphere
		len = size;
		if ((int)cell->endN.z<1234 && vdif(cell->endN,cell->endC)>NOISE) { // poles are set
			if (cell->child[0]->cent.x < 0.0) {	// but kids have no cent.[xyz] distances
				FOR(i,n) { Cells *kidi = cell->child[i]; // set cent reference distances
					kidi->cent.x = vdif(kidi->xyz,cell->endN);
					kidi->cent.y = vdif(kidi->xyz,cell->xyz);
					kidi->cent.z = vdif(kidi->xyz,cell->endC);
				}	
			}
		} else { return; } // poles are not set
	}
	if (type==2)	// tube
	{ float fn = (float)n;
		if (moltype==0) { // protein
			if (sort==0) len = size*sqrt(fn)*loopAXIS;	// loops
			if (sort==1) len = size*fn*bondCA*alphAXIS/3.8;	// alpha
			if (sort==2) len = size*fn*bondCA*betaAXIS/3.8;	// beta
		}
		if (moltype==1) { // nucleic
			if (sort==0) len = 0.2*sqrt(fn);	// loop
			// RNA = 2.3 rise/bp --> 0.115 (1/2 for duplex, 1/10 for scale)
			if (sort==1) len = 0.12*fn;		// stem
		}
	}
	if (type==3)	// ellipsoid
	{ float fn = E/2, nth = 1.0/fn, a =log2(fn+1),
	  	f = (fn-1.0)/(pow(fn,a)*(pow(2.0,a)-1.0));
       		len = size*(nth+f*pow((float)sort,a)); // set axis length for ellipsoid
	}
	len *= 0.5; // used as semiaxis length
	set = 0;
	if (cell->child[0]->cent.x > 0.0) set = 1; // centre distances exist
	if (n == 1) { // simple fix for tube/ellip with just one child
		if (type > 1) vave(cell->endN,cell->endC, &(cell->xyz)); // spheres keep their centre
		set = -1;
	}
	if (level<depth-1 && set > 0) { // reset axes (using cent distances)
		setWcent(cell); // set cell position to centroid of weighted kids
		inN = inC = 0;
		vinit(&newN); vinit(&newC);
		FOR(i,n)
		{ Cells *kidi = cell->child[i];
		  float a = kidi->cent.x, b = kidi->cent.y, c = kidi->cent.z;
			if (goodTri(a,b,len)) {
				findPole(cell->xyz,cell->endN,kidi->xyz,a,b,len,&new);
				vsum(new, &newN);
				inN++;
			}
			if (goodTri(b,c,len)) {
				findPole(cell->xyz,cell->endC,kidi->xyz,c,b,len,&new);
				vsum(new, &newC);
				inC++;
			}
		}
		set = 0;
		if (inN && inC)  {
			// average the new and old positions 
			vdiv(&newN,(float)inN); vsub(newN,cell->endN, &shift);
			vmul(&shift,weight); vadd(cell->endN,shift, &(cell->endN));
			vdiv(&newC,(float)inC); vsub(newC,cell->endC, &shift);
			vmul(&shift,weight); vadd(cell->endC,shift, &(cell->endC));
			if (type > 1) vave(cell->endN,cell->endC, &(cell->xyz)); // set tube and ellipsoid midpoint 
			set = 1;
		}
		if (set) {
			// reset N+C separation to exact length along the new axis
			vatA(cell->xyz,cell->endN, &(cell->endN),len);
			vatA(cell->xyz,cell->endC, &(cell->endC),len);
		}
	}
	FOR(i,n) fixANYaxis(cell->child[i], level+1, weight);
}

void findPole ( Vec vA, Vec vE, Vec vB, float a, float b, float g, Vec *new )
{
float	h, r;
Vec	vC,vD,vF,vG, vAB,vFD;
//	a = kidi->cent.[xz] = start d to pole (end[NC])
//	b = kidi->cent.y    =  start d to centre (xyz)
	if (b<NOISE) { // B is on A
		vcopy(vE,new);
		return;
	}
//
//                  G<-E
//                 .C  : 
//            g .  /|  :
//           .   a/ |r :
//        .      /  |  :
//      A-------B---D..F     
//          b     h  
//
// g=len/2 (half the pole-pole axis length)
// A = cell->xyz, B = kidi->xyz, C = pole position to be estimated.
// 1) find the radius r of the circle at B normal to AB such that:
// r2 = a2 - h2 = g2 - (b+h)2
//              = g2 - b2 -2bh - h2
//          2bh = g2 - b2 - a2
//            h = (g2-b2-a2)/(2b) 
//            r = sqrt(a2-h2)
// 2) find the centre of the circle D distance h beyond B
// 3) find the closesest point on the circle to the old pole, E = cell->end[NC]
// 	a) get F = image of E on the circle axis AD
// 	b) shift E to G (by FD)
// 	c) find C at distance r along DG
//
	h = (g*g-b*b-a*a)/(2.0*b);
	r = sqrt(a*a-h*h);
	vatA(vA,vB,&vD,b+h);
	dot2line(vA,vD,vE,&vF);
	vsub(vD,vF,&vFD);
	vadd(vE,vFD, &vG);
	vatA(vD,vG,&vC,r);
	vcopy(vC,new);
}	

randomise( Cells *cell, int level )
{
int	i;
	///cell->xyz.z *= -1.0;
	///vradd(&(cell->xyz),5.0);
	///if (level==2) vradd(&(cell->xyz),5.0);
	///if (level==2) vrset(&(cell->xyz),1.0);
	///vrset(&(cell->xyz),10.0);
	for (i=0; i<cell->kids; i++) randomise(cell->child[i], level+1);
}
///////////
for (i=0; i<cell->kids; i++) { Cell *ci = cell->child[i];
printf("ATOM%7d  CA  GLY A%4d     %7.3f %7.3f %7.3f  0.00 %5.2f\n", i,i,ci->xyz.x,ci->xyz.y,ci->xyz.z,(float)ci->lost);
} i=0;
printf("TER ATOM\n");
printf("ATOM%7d  CA  CYS A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,cell->endN.x,cell->endN.y,cell->endN.z);
printf("ATOM%7d  CA  CYS A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,cell->endC.x,cell->endC.y,cell->endC.z);
printf("TER ATOM\n");
exit(1);
///////////
*/
