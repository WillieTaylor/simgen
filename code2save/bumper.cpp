#include "util.hpp"
#include "geom.hpp"
#include "cell.hpp"
#include "data.hpp"

float inEgg ( Vec, Seg, float );

void bumper ()
{
	Cell::world->bumps();
}

bool exempt ( Cell *a, Cell *b )
{ // returns 1 if <a> and <b> are exempt from bumping
//if (a->hit) { Pi(a->level) Pi(a->id) Pi(a->bump) NL }
/*
	if (a->hit && a->hit==b) return 1;
	if (b->hit && b->hit==a) return 1;
*/
	FOR(k,a->nbonds) { if (a->bond[k].to==b) return 1; }
	FOR(k,b->nbonds) { if (b->bond[k].to==a) return 1; }
//	FOR(k,a->nlinks) { if (a->link[k].to==b) return 1; }
//	FOR(k,b->nlinks) { if (b->link[k].to==a) return 1; }
	return 0;
}

#define SWITCH 20	// number of cells below which bump testing is pairwise
#define MAXHOLD 9999

int sortBumps ( const void *ac, const void *bc )
{
        Bumps   *a = (Bumps*)ac, *b = (Bumps*)bc;
        if (a->d < b->d) return -1;
        if (a->d > b->d) return  1;
        return 0;
}

int getBumpin ( Cell *cell, Bumps *pairs ) {
// return the number (and <list>) of possible bumping pairs in cell <c> closer than <bump>
// <bump>: spheres=rad.s, tubes=lengths, ellipsoids=max-axes
Cell	*a, *b;
float	alen, blen;
float	d, dd, bump, bbump;
Cell	*kid = cell->child[0];
Data    *p = Data::model+kid->model;
int	moltype = p->moltype;
int	kidlev = kid->level;
float	size = p->bumps[kidlev];
int	type = cell->type,
	kids = cell->kids,
	hold = HOLD*kids,
	i, j, k, m, n=0,
	flat = 0;
int	close[MAXHOLD];
	if (hold > MAXHOLD) hold = MAXHOLD;
	if (kids==1) return 0;
	if (moltype==3 && Data::shrink > NOISE) flat = 1; // avoid Z ranks for 2D cells
	if (size<0) size = -size;
	if (kids < SWITCH || flat) {	// use pairwise
		for (i=0; i<kids-1; i++) {
			a = cell->child[i];
			if (a->empty) continue;
			for (j=i+1; j<kids; j++) {
				b = cell->child[j];
				if (b->empty) continue;
				if (kidlev==Data::depth && exempt(a,b)) continue; // just filter atom level
				alen = a->len; blen = b->len;
				if (a->type==2) alen += size;	// add size of end-caps
				if (b->type==2) blen += size;	// to both tube objects
				if (alen > size) bump = alen; else bump = size;
				if (blen > size) bump += blen; else bump += size;
				bump *= 0.5; bbump = bump*bump;
				dd = a->xyz||b->xyz;
				if (dd > bbump) continue;
				pairs[n].c = 0; // not used
				pairs[n].a = a; pairs[n].b = b;
				pairs[n].d = sqrt(dd);
				pairs[n].bump = bump;
				n++;
				if (n == hold) break;
			}
			if (n == hold) break;
		}
		if (n==0) return 0;
		qsort(pairs,n,sizeof(Bumps),sortBumps);         // sort on separation (closest first)
		return n;
	} else	// use sorted lists to find pairs
	{ int   got, span = 5+(int)sqrt((float)kids);
		for (i=0; i<kids; i++) {
			a = cell->child[i];
			if (a->empty) continue;
			m = 0;
			for (j=0; j<3; j++) // gather locally ranked children in each dimension
			{ int	rat = a->ranks[j], out;
				out = 0;
				for (k=rat-span; k<rat+span; k++) {	// check the local rank list
					if (k>=kids) break;		// intra-family so max=kids
					if (k<0) continue;		// not on list yet
					if (k==rat) continue;		// skip self
					if (n == hold) break;
					b  = cell->rank[k][j];
					if (b->empty) continue;
					if (kidlev==Data::depth && exempt(a,b)) continue; // filter atom level
					close[m] = 10*(b->uid) + j + 1;	// code dim as j/10
					m++;
					if (m==hold) break;
				}
				if (m==hold) break;
			}
			sort(close,m);
			got = 1; // 123 = XYZ (in sorted order in pairs)
			FOR(ii,m) { int hit, dim;  // keep only children present in all 3 dimensions
				hit = (int)(0.1*(float)close[ii]);
				dim = close[ii]-hit*10;
				if (dim != got) continue;
				if (got==3) { float d; // keep
					got = 0;
					b = Cell::uid2cell[hit];
					if (a->len > size) bump = a->len; else bump = size;
					if (b->len > size) bump += b->len; else bump += size;
					bump *= 0.5; bbump = bump*bump;
					dd = a->xyz||b->xyz;
					if (dd < bbump) {
						pairs[n].c = 0; // not used
						pairs[n].a = a; pairs[n].b = b;
						pairs[n].d = sqrt(dd);
						pairs[n].bump = bump;
						n++;
					}
					if (n==hold) break;
				}
				if (n==hold) break;
				got++;
			}
		}
		if (n==0) return 0;
		qsort(pairs,n,sizeof(Bumps),sortBumps);         // sort on separation (closest first)
		m = 1;
		FOR(ii,n) { // remove duplicate entries
			if (ii==0) continue;
			if ((pairs[ii].a->uid==pairs[ii-1].a->uid) && (pairs[ii].b->uid==pairs[ii-1].b->uid)) continue;
			if ((pairs[ii].a->uid==pairs[ii-1].b->uid) && (pairs[ii].b->uid==pairs[ii-1].a->uid)) continue;
			pairs[m] = pairs[ii];
			m++;
		}
		return m;
	}
}

float tube_to_egg ( Cell *a, Cell *b ) {
// returns an approximation to the closest approach of a tube <a> to an ellipsoid <b> surface
// NB the value returned by inEgg() is not a surface distance but is zero on the surface
// NB assumes radially symmeric ellipsoid
Seg	ca = Seg(a->endN,a->endC),	cb = Seg(b->endN,b->endC);
Data	*pa = Data::model+a->model,	*pb = Data::model+b->model;
float	sizea = pa->sizes[a->level],	sizeb = pb->sizes[b->level];
float	d1,d2,d3;
Vec	p1,p2,p3;
	if (a->type!=2 && b->type!=3) { Pt(Bad types in tube_to_egg()) Pi(a->type) Pi(b->type) NL exit(1); }
	p1 = a->endN; p2 = a->xyz; p3 = a->endC;
	DO { float d;
		d1 = inEgg(p1,cb,sizeb);
		d2 = inEgg(p2,cb,sizeb);
		d3 = inEgg(p3,cb,sizeb);
		if (d1<0 || d2<0 || d3<0) return -999.9; // flag bump;
		d = p1|p3;
		if (d < 0.01) { // pretty close
			d = vec_to_egg(p2,cb,sizeb);
			return d - sizea*0.5;
		}
		if (d1+d2 < d2+d3) {
			p3 = p2; p2 = p1 & p2;
		} else {
			p1 = p2; p2 = p3 & p2;
		}
	}
}

// tested in newsims/bumpell
#define Nrand  10 // 20 // WAS 500 for testing
#define Nhold 500

Seg trisect ( int, Vec*, Vec*, Vec*, Vec*, Vec, Vec, Mat, Mat, Vec, Vec, Vec, Vec, float, int );
Vec shell ( Vec, Vec, Vec, Mat, Vec );
Vec shall ( Vec, Vec, Vec, Mat, Vec );
Vec sholl ( Vec, Vec, float, float, Vec );
float bumping ( int, Vec*, Vec*, Vec, Vec, Vec, Vec, Mat, Mat );
float in_egg ( Vec, Vec, Vec, Mat );

float egg_to_egg ( Cell *a, Cell *b ) {
// returns a good approximation to the closest approach of two (non-scalene) ellipsoid surfaces
Data	*pa = Data::model+a->model,	*pb = Data::model+b->model;
float	sizea = pa->sizes[a->level],	sizeb = pb->sizes[b->level],	// A=B diameter
	denda = a->endN|a->endC,	dendb = b->endN|b->endC,	// C axis length
	dmina,dminb, dmaxa,dmaxb, dmin,dmax, da,db, d;
Seg	ca = Seg(a->endN,a->endC), cb = Seg(b->endN,b->endC), c, p;
Seg	ends, best;
float	dends, wsum;
float	**mat, score, rand = 0.2; // WAS 0.5 for testing
Mat	axesA, axesB;
Vec	axisA, axisB, centA, centB, bestA, bestB, surfA, surfB;
Vec	quad[8], quadA[8], quadB[8], sectA[6], sectB[6];
Vec	dispA[Nrand], dispB[Nrand];
Vec	mid, ave, aimA, aimB;
Vec	moveA, moveB;
Pairs	dist[64];
int	rank[64];
int	qA,qB, qAtop,qBtop, n = 0;
Vec	allA[Nhold], allB[Nhold];
Vec	randout;
int	on = 0;
	if (a->ends < 0) sizea = denda/Data::Eratio[a->sort]; // use given ratio and C length
	if (b->ends < 0) sizeb = dendb/Data::Eratio[b->sort]; // use given ratio and C length
	dmina = min(sizea,denda);	dminb = min(sizeb,dendb);
	dmaxa = max(sizea,denda);	dmaxb = max(sizeb,dendb);
	dmin = (dmina+dminb)*0.5;	dmax = (dmaxa+dmaxb)*0.5;
	d = a->xyz|b->xyz;
	if (d > dmax) return 999.9;	// beyond maximum contact distance
	if (d < dmin) return d-dmin;	// centres closer than minimum separation
	sizea *= 0.5; sizeb *= 0.5; denda *= 0.5; dendb *= 0.5;	// reset to semiaxis lengths
	mat = new float*[6]; FOR(i,6) mat[i] = new float[6];
	centA = a->xyz; centB = b->xyz;
	mid = (centA+centB)*0.5;
	axisA.x = denda;  axisA.y = axisA.z = sizea;
	axisB.x = dendb;  axisB.y = axisB.z = sizeb;
	axesA.A = a->endC - a->xyz;
	randout = (axesA.A^centB).norm();
	axesA.B = (axesA.A^randout).getVec(sizea);
	axesA.C = (axesA.A^axesA.B).getVec(sizea);
	axesB.A = b->endC - b->xyz;
	axesB.B = (axesB.A^randout).getVec(sizeb);
	axesB.C = (axesB.A^axesB.B).getVec(sizeb);
	quad[0] =  Vec( 1, 1, 1); quad[1] =  Vec( 1, 1,-1); 
	quad[2] =  Vec( 1,-1, 1); quad[3] =  Vec( 1,-1,-1);
	quad[4] =  Vec(-1, 1, 1); quad[5] =  Vec(-1, 1,-1);
	quad[6] =  Vec(-1,-1, 1); quad[7] =  Vec(-1,-1,-1);
	FOR(i,8) // make real axis end-points for both ellipsoids and put surface centroids in quad[AB]
	{ Vec	midA, midB;
		midA = centA + (axesA.A*quad[i].x + axesA.B*quad[i].y + axesA.C*quad[i].z)/3.0;
		allA[i] = quadA[i] =  shell(midA,centA,axisA,axesA,quad[i]);
		midB = centB + (axesB.A*quad[i].x + axesB.B*quad[i].y + axesB.C*quad[i].z)/3.0;
		allB[i] = quadB[i] =  shell(midB,centB,axisB,axesB,quad[i]);
	}
	sectA[0] = centA+axesA.A; sectB[0] = centB+axesB.A;
	sectA[1] = centA+axesA.B; sectB[1] = centB+axesB.B;
	sectA[2] = centA+axesA.C; sectB[2] = centB+axesB.C;
	sectA[3] = centA-axesA.A; sectB[3] = centB-axesB.A;
	sectA[4] = centA-axesA.B; sectB[4] = centB-axesB.B;
	sectA[5] = centA-axesA.C; sectB[5] = centB-axesB.C;
	FOR(i,6) { allA[i+8] = sectA[i]; allB[i+8] = sectB[i]; }
	d = bumping(14, allA,allB, centA,centB, axisA,axisB, axesA,axesB);
	if (d>NOISE) return -d;
	n = 0;
	dends = 9999.9;
	wsum = 0.0;
	ave.zero();
	d = centA | centB;
	FOR(i,14) FOR(j,14)
	{ float w, dij; Vec aij; Vec pAi = allA[i], pBj = allB[j];
		dij = pAi|pBj;  aij = pAi+pBj;
		if (dij > d) continue;
		if (dij < dends) { dends = dij; best.A = pAi; best.B = pBj; }
		w = 1.0/(dij*dij); wsum += w;
		ave += aij*0.5*w;
		n++;
	}
	ave /= wsum;
	aimA = shell(ave,centA,axisA,axesA,quad[0]);
	aimB = shell(ave,centB,axisB,axesB,quad[0]);
	n = 0;
	FOR(i,8) FOR(j,8) {
		dist[n].a = i; dist[n].b = j; 
		dist[n].s = (quadA[i]|quadB[j])+(aimA|quadA[i])+(aimB|quadB[j]);
		n++;
	}
	sort(dist,rank,-n);
	dmin = 999.9;
	allA[6] = aimA; allB[6] = aimB;
	FOR(i,3) // test just the top 3
	{ int	r = rank[i], qA = dist[r].a, qB = dist[r].b;
	  float dab;
		// pass quadrant end-points and mid-points in sect[AB] to trisect()
		sectA[0] = centA + axesA.A*quad[qA].x;
		sectA[2] = centA + axesA.B*quad[qA].y;
		sectA[4] = centA + axesA.C*quad[qA].z;
		sectA[1] = shell((sectA[0]+sectA[2])*0.5,centA,axisA,axesA,quad[qA]);
		sectA[3] = shell((sectA[2]+sectA[4])*0.5,centA,axisA,axesA,quad[qA]);
		sectA[5] = shell((sectA[4]+sectA[0])*0.5,centA,axisA,axesA,quad[qA]);
		sectB[0] = centB + axesB.A*quad[qB].x;
		sectB[2] = centB + axesB.B*quad[qB].y;
		sectB[4] = centB + axesB.C*quad[qB].z;
		sectB[1] = shell((sectB[0]+sectB[2])*0.5,centB,axisB,axesB,quad[qB]);
		sectB[3] = shell((sectB[2]+sectB[4])*0.5,centB,axisB,axesB,quad[qB]);
		sectB[5] = shell((sectB[4]+sectB[0])*0.5,centB,axisB,axesB,quad[qB]);
		FOR(j,6) { allA[j] = sectA[j]; allB[j] = sectB[j]; }
		ends = trisect(7,allA,allB,sectA,sectB,centA,centB,axesA,axesB,axisA,axisB,quad[qA],quad[qB],9999.9,0);
		if (ends.B.z > 9999.0) return -ends.B.x;	// found bumping
		dab = ends.len();
		if (dab > dmin) continue;
		best = ends; dmin = dab; qAtop = qA; qBtop = qB; n = i;
	}
	rand *= sqrt(dmin);
	moveA = best.A;
	FOR(i,Nrand) { Vec v = moveA; v.set_disp(rand);
		surfA = shell(v,centA,axisA,axesA,quad[qAtop]);
		dispA[i] = surfA;
		d = surfA|best.B; if (d<dmin) { moveA = surfA; dmin = d; }
	}
	moveB = best.B;
	FOR(i,Nrand) { Vec v = moveB; v.set_disp(rand);
		surfB = shell(v,centB,axisB,axesB,quad[qBtop]);
		dispB[i] = surfB;
		d = surfB|moveA; if (d<dmin) { moveB = surfB; dmin = d; }
	}
	FOR(i,Nrand) FOR(j,Nrand) {
		d = dispA[i]|dispB[j];
		if (d<dmin) {
			moveA = dispA[i]; moveB = dispB[j]; dmin = d;
		}
	}
	da = vec_to_egg(moveA,b->endN,b->endC,sizeb*2.0); if (da<dmin) dmin = da;
	db = vec_to_egg(moveB,a->endN,a->endC,sizea*2.0); if (db<dmin) dmin = db;
	if (da>0 && db>0) { // monitor accuracy and print errors over 50%
		d = fabs(da-db)/(1.0+da+db);
		if (d>0.5) { Pt(Poor egg2egg estimate) Pr(d) Pr(da) Pr(db) NL }
	}
	return dmin;
}

float check_normal ( Vec A, Vec a0, Vec a1, Vec a2,  Vec B, Vec b0, Vec b1, Vec b2 ) {
// return the sum of the distances that the normals from <A> and <B> pass their opposing points
Vec	normA, normB;
float	toA, toB;
	normA = A+((a1-a0)^(a2-a0));
	normB = B+((b1-b0)^(b2-b0));
	toA = A.vec_to_line(B,normB);
	toB = B.vec_to_line(A,normA);
	return toA+toB;
}

Seg trisect ( int on, Vec *allA, Vec *allB, Vec *sectA, Vec *sectB, Vec centA, Vec centB, Mat axesA, Mat axesB, Vec axisA, Vec axisB, Vec qA, Vec qB, float last , int depth ) {
int	set[4][3] = { 1,2,3, 3,4,5, 5,0,1, 1,3,5 };
Vec	midA,midB, bestA,bestB, surfA,surfB;
int	triA0,triB0, triA1,triB1, triA2,triB2;
float	wsum, w, d, dif, dab, dmin;
Vec	A0,A1,A2, B0,B1,B2;
Vec	aim, aimA, aimB;
Vec	a[8], b[8];
Seg	ends;
int	p,q;
	aimA.zero();
	aimB.zero();
	wsum = 0.0;
	FOR(i,on) FOR(j,on) {
		d = allA[i]|allB[j]; w = 1.0/(d*d);
		aimA += allA[i]*w; aimB += allB[j]*w;
		wsum += w;
	} // set target points <aimA>,<aimB> to weighted mean and project onto surface
	aimA /= wsum; aimB /= wsum;
	aimA = shell(aimA,centA,axisA,axesA,qA);
	aimB = shell(aimB,centB,axisB,axesB,qB);
	allA[on] = aimA; allB[on] = aimB; on++;	// append to current surface points
	dmin = 9999.9;
	FOR(i,on) FOR(j,on) { // find min of all current surface pairs
		d = allA[i]|allB[j];
		if (d<dmin) { dmin = d; p = i; q = j; }
	}
	aimA = allA[p]; aimB = allB[q];
	FOR(i,4) // make midpoints for the 4 triangles in segments A and B
	{ int	si0 = set[i][0], si1 = set[i][1], si2 = set[i][2];
		midA = (sectA[si0]+sectA[si1]+sectA[si2])/3.0;
		allA[on+i] = a[i] = shell(midA,centA,axisA,axesA,qA);
		midB = (sectB[si0]+sectB[si1]+sectB[si2])/3.0;
		allB[on+i] = b[i] = shell(midB,centB,axisB,axesB,qB);
	}	// and append to collection of points <on> the surfaces 
	on += 4;
	d = bumping(on, allA,allB, centA,centB, axisA,axisB, axesA,axesB);
	if (d>NOISE) { 
		return Seg(Vec(d,0,9999.9)); // Seg.B.z > 9999 = bump 
	}
	dmin = 9999.9;
	w = (float)depth;
	if (depth > 1) w = sqrt(w);
	FOR(i,4) // loop over the 4 midpoints in segment A
	{ int	si0 = set[i][0], si1 = set[i][1], si2 = set[i][2];
		midA = a[i];
		FOR(j,4) // loop over the 4 midpoints in segment B
		{ int	sj0 = set[j][0], sj1 = set[j][1], sj2 = set[j][2];
			midB = b[j];
			d = check_normal(midA,sectA[si0],sectA[si1],sectA[si2],midB,sectB[sj0],sectB[sj1],sectB[sj2]);
			dab = w*sqrt(d)+(2.0/(1+w))*(midA|midB)+(aimA|midA)+(aimB|midB);
			if (dab<dmin) {
				dmin = dab;
				triA0=si0; triA1=si1; triA2=si2;
				triB0=sj0; triB1=sj1; triB2=sj2;
				p = i; q = j;
			}
		}
	}
	// set best triangle pair for next iteration
	A0 = sectA[triA0]; A1 = sectA[triA1]; A2 = sectA[triA2];
	B0 = sectB[triB0]; B1 = sectB[triB1]; B2 = sectB[triB2];
	surfA = a[p]; surfB = b[q];
	dab = surfA|surfB;
	dif = last-dab;
	//if ((depth > 3 && last-dab < 0.001)||(depth > 6)) {
	if (depth==3) { // go with whatever (4xfaster than above)
		// NB for depth > 7, grid points get too close. For depth=9, on=92 (check Nhold)
		d = bumping(on, allA,allB, centA,centB, axisA,axisB, axesA,axesB);
		if (d>NOISE) { 
			return Seg(Vec(d,0,9999.9)); // Seg.B.z > 9999 = bump 
		}
		dmin = 9999.9;
		FOR(i,on) FOR(j,on) { // find min of all current surface pairs
			d = allA[i]|allB[j];
			if (d<dmin) { dmin = d; p = i; q = j; }
		}
		surfA = allA[p]; surfB = allB[q];
		return Seg(surfA,surfB);
	}
	sectA[0] = A0; allA[on+0] = sectA[1] = shell((A0+A1)*0.5,centA,axisA,axesA,qA);
	sectA[2] = A1; allA[on+1] = sectA[3] = shell((A1+A2)*0.5,centA,axisA,axesA,qA);
	sectA[4] = A2; allA[on+2] = sectA[5] = shell((A2+A0)*0.5,centA,axisA,axesA,qA);
	sectB[0] = B0; allB[on+0] = sectB[1] = shell((B0+B1)*0.5,centB,axisB,axesB,qB);
	sectB[2] = B1; allB[on+1] = sectB[3] = shell((B1+B2)*0.5,centB,axisB,axesB,qB);
	sectB[4] = B2; allB[on+2] = sectB[5] = shell((B2+B0)*0.5,centB,axisB,axesB,qB);
	d = bumping(on+3, allA,allB, centA,centB, axisA,axisB, axesA,axesB);
	if (d>NOISE) { 
		return Seg(Vec(d,0,9999.9)); // Seg.B.z > 9999 = bump 
	}
	// pass sub-triangles (and append to surface collection)
	ends = trisect(on+3,allA,allB,sectA,sectB,centA,centB,axesA,axesB,axisA,axisB,qA,qB,dab,depth+1);
	return ends;
}

float bumping ( int n, Vec *allA, Vec *allB, Vec centA, Vec centB, Vec axisA, Vec axisB, Mat axesA, Mat axesB )
{
int	p, q = 0;
float	d = 0.0, dmin = 999.9;
Vec	surfA, surfB, quad = Vec(1,1,1);
	FOR(i,n) // check if any points lie inside the other ellipsoid
	{ float din;
		din = in_egg(allA[i],centB,axisB,axesB);
		if (din < 1.0-NOISE) { // A is in B (flag with q=1)
			if (din<dmin) { dmin = din; p = i; q =  1; }
		}
		din = in_egg(allB[i],centA,axisA,axesA);
		if (din < 1.0-NOISE) { // B is in A (flag with q=-1)
			if (din<dmin) { dmin = din; p = i; q = -1; }
		}
	}
	if (q>0) { // A is deeper inside B
		surfB = shell(allA[p],centB,axisB,axesB,quad);
		d = surfB|allA[p];
		// Pi(n) Pr(dmin) Pt(A inside B by) Pr(d) NL
	}
	if (q<0) { // B is deeper inside A
		surfA = shell(allB[p],centA,axisA,axesA,quad);
		d = surfA|allB[p];
		// Pi(n) Pr(dmin) Pt(B inside A by) Pr(d) NL
	}
	return d;
}

float in_egg ( Vec p, Vec c, Vec x, Mat axes ) {
// Returns the fractional distance <d> of point <p> to the ellipsoid surface.
// <d> < 1 = inside, <d> > 1 = outside
// The ellipsoid has axes <M> (at 0 with lengths x)
float	d;
Vec	r, q = (p-c).get_frac(axes);		// q = components of <p> in basis set <M>
	r = Vec(q.x*x.x, q.y*x.y, q.z*x.z);	// r = real-space position of <p> in <M>
	d = r.x*r.x/(x.x*x.x) + r.y*r.y/(x.y*x.y) + r.z*r.z/(x.z*x.z); // surface definition
	return d;
}

Vec shall ( Vec line, Vec cent, Vec axis, Mat axes, Vec quad ) {
// returns the point on the ellipsoid (<axes>=ABC at 0) surface cut by a <line> from the <cent>re
/*
if the axes lie on XYZ then,
	xx/AA + yy/BB + zz/CC = 1
	x = pa, y = pb, z = pc
	pp = 1/(aa/AA+bb/BB+cc/CC)
	line = abc, surf = xyz
*/
Mat	frame;
Vec	surf;
float	A,B,C, a,b,c, p;
	line -= cent;					// quadrant centroid shifted to origin
	frame = Mat(axes.A.getVec()*quad.x, axes.B.getVec()*quad.y, axes.C.getVec()*quad.z);
	surf = line*frame;				// rotated to A=X, B=Y, C=Z
	A=axis.x, B=axis.y, C=axis.z;			// semi-axis lengths
	a=surf.x, b=surf.y, c=surf.z;			// rotated point
	p = 1.0/sqrt(a*a/(A*A)+b*b/(B*B)+c*c/(C*C));	// scale factor
	surf *= p;					// to put point on ellipsoid surface
	surf *= frame.get_trans();			// rotated back to ABC frame
	return cent+surf;				// added back to centre
}

Vec sholl ( Vec line, Vec cent, float A, float B, Vec axis ) {
// returns the point on the ellipsoid (<axes> = A>B=C at 0) surface cut by a <line> from the <cent>re
/*
in the plane of the major axis (A) and the <line> with components a,b to A,
the point where the line cuts the surface has corresponding components g,h.
Now	gg/AA + hh/BB = 1
and	g/a = h/b
so	gg = AA(1-hh/BB) = aa.hh/bb
	AA - hh.AA/BB = hh.aa/bb
	AA = hh.aa/bb + hh.AA/BB
	hh = AA/(aa/bb+AA/BB)
*/
Vec	surf;
float	AA=A*A, BB=B*B, aa,bb,b, d,ff,gg,hh;
	line -= cent;					// shift line to origin
	b = line.vec_to_line(axis);			// perpendicular dist from line to axis
	bb = b*b;
	ff = line.sqr();
	aa = ff-bb;
	hh = AA/(aa/bb+AA/BB);
	gg = aa*hh/bb;
	d = sqrt((gg+hh)/ff);
	surf = line*d;					// extend <line> to ellipsoid surface
	return cent+surf;				// added back to centre
}

Vec shell ( Vec line, Vec cent, Vec axis, Mat axes, Vec quad ) {
float	ab, bc, ca, close = NOISE;
	ab = axis.x - axis.y;
	bc = axis.y - axis.z;
	ca = axis.z - axis.x;
	// check for oblate/prolate
	if (ab*ab < close) return sholl(line,cent,axis.z,axis.x,axes.C);
	if (bc*bc < close) return sholl(line,cent,axis.x,axis.y,axes.A);
	if (ca*ca < close) return sholl(line,cent,axis.y,axis.z,axes.B);
	// otherwise scalene
	return shall(line,cent,axis,axes,quad);
}
