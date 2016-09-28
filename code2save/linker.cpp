#include "util.hpp"
#include "geom.hpp"
#include "cell.hpp"
#include "data.hpp"

#define BOND 1.0 // 0.3
#define LINK 1.0 // 0.6

void bondStem ( Cell* );
void fixSheet ( Cell* );

Vec	rdisp ( Vec, Vec, float );
float	scoreFit ( Vec, Vec, Vec, Vec, Vec, Vec, float );
void	fixHinge ( Cell*, Bonds*, float );

void bonder ( Cell *cell )
{
float	d;
Vec	shift;
Cell	*a, *b;
int	i, j, n, m = cell->model;
int	moltype, subtype, level = cell->level;
float	sislen, brolen, size, basic, kick = BOND;
float	*sizes = Data::model[m].sizes;
float	*bonds = Data::model[m].bonds;
int	*chain = Data::model[m].chain;
int	*local = Data::model[m].local;
float	*bondlen;
char	**bondstr;
int	nbonds;
int	end = 0;
	if (cell==0) return;
	if (cell->empty) return;
	n = cell->kids;
	if (level && cell->parent->solid > 0) return; 
	moltype = Data::model[m].moltype;
	subtype = Data::model[m].subtype;
	size = sizes[level]+bonds[level];
	if (level && cell->sis && (cell->sis->level != cell->level)) end = -1; // N end of a chain
	if (level && cell->bro && (cell->bro->level != cell->level)) end =  1; // C end of a chain
	if (level && chain[level]) { // in a chain
		sislen = brolen = size;
		if (moltype==0) {  // PROTein chain (not on an end)
			if (local[level]) {  // using preset lengths (set by -ve size in param.model
				if (end>-1) sislen = cell->prox.x;	// use specific sis length (x)
				if (end< 1) brolen = cell->prox.z;	// use specific bro length (z)
			}
			if (sislen >  999.0) sislen = size;	// use generic bond length
			if (brolen >  999.0) brolen = size;	// dito for bro
			if (level < depth) {
				// check if bonded at lower level (-ve sis/brolen = skip)
				if (cell->starts->sis->level != level) sislen = -1.0;
				if (cell->finish->bro->level != level) brolen = -1.0;
                		// use HINGE data for protein or nucleic
                        	FOR(i,cell->nbonds) { Bonds *hinge = cell->bond+i;
                                	if (hinge->to == 0) continue; // no bond
					if (hinge->link==0) continue; // normal bond
                                	fixHinge(cell,hinge,bonds[level]);
					sislen = brolen = -1.0; // dont' refine centroids
                        	}
			}
			if (level==depth-1) { // dont bond loop tubes with generic length
				if ((cell->sort==0 || cell->sis->sort==0) && sislen==size) sislen = -1.0;
				if ((cell->sort==0 || cell->bro->sort==0) && brolen==size) brolen = -1.0;
			}
		}
		if (moltype==1) { // Nucleic acid chain
			if (level==depth-1) {	// don't refine basepaired stack bond distance in a loop
				sislen = brolen = 0.5;
				if (cell->sort==0) sislen = brolen = -1.0;
				if (cell->sis->sort==0) sislen = -1.0;
				if (cell->bro->sort==0) brolen = -1.0;
			}
			if (level==depth-2 && cell->type==2) {	// bond domain tube ends with 0 length bond
				sislen = brolen = -1.0;		// skips call of part2cells() below
				if (subtype==1) { Vec v; float w = 0.1;	// link DNA (but don't join RNA ends)
					if (cell->sis->level==level) { Vec Nter, Cter; // fix sis side
						v = (cell->sis->endC - cell->sis->endN).norm();
						Cter = cell->sis->xyz + v*(cell->sis->len*0.6);
						v = (cell->endN - cell->endC).norm();
						Nter = cell->xyz + v*(cell->len*0.6);
						v = (Nter-Cter)*w;
						cell->move(-v); cell->sis->move(v);
					}
				}
				if (subtype==0) // bond RNA stem ends
				{ Cell	 *wC, *cC, *cN, *wN=0;
					// find stem end points
					FOR(i,cell->kids) { Cell *kidi = cell->child[i];
						if (kidi->sort==0) continue;
						if (wN==0) {
							wN = kidi->child[0];
							cN = kidi->child[1];
						}
						wC = kidi->child[0];
						cC = kidi->child[1];
					}
					// refine end distances
					bondStem(wN);
					bondStem(cN);
					bondStem(wC);
					bondStem(cC);
					sislen = brolen = -1.0;  // no more refinement
				}
			}
		}
		if (cell->bond[0].link < -990) sislen = -1.0;
		if (cell->bond[1].link < -990) brolen = -1.0;
		// sis/brolen<0 = flagged by RELINK as not for refinement
		if (randf()<0.5) {	// fix sis side
			if (end>-1 && sislen>0.0) part2cells(cell,cell->sis,sislen,kick);
		} else {		// fix bro side
			if (end< 1 && brolen>0.0) part2cells(cell,cell->bro,brolen,kick);
		}
	}
	nbonds = abs(chain[level]);
	if (nbonds > 1) { // may be cross-linked
		for (int i=2; i<=nbonds; i++) {
			if (cell->bond==0 || cell->bond[i].to==0) continue;
			part2cells(cell,cell->bond[i].to,sislen,kick);
		}
	}
	if (randf()<0.5) { // randomly alternate direction
		FOR(i,n) bonder(cell->child[i]);
	} else {
		FOR(i,n) bonder(cell->child[n-i-1]);
	}
}

void bondStem (Cell *at) {
Vec	v;
float	d, e, bond = 3.5;
int	atend, toend, ends;
Cell	*ats = at->parent->parent, // stem
	*to, *tob, *tos;
	// check if at links (via loop) to a different stem 
	to = at->bro;
	if (to->level != depth) return;		// end of chain
	DO { // follow loop to its end (loop has sort=0, basepair and stem levels have sort=1)
		tos = to->parent->parent;
		if (to->parent->sort == 1) break; // in another basepair
		if (tos->uid != ats->uid) break;  // in a different stem
		to = to->bro;
	}
	if (tos->uid == ats->uid) return;	// still in same stem
	if (tos->level != ats->level) return;	// shouldn't happen
	// Pi(at->resn) Pi(to->resn) Pi(ats->uid) Pi(tos->uid) 
	if ((at->xyz|ats->endN) < (at->xyz|ats->endC)) atend = 0; else atend = 1;
	if ((to->xyz|tos->endN) < (to->xyz|tos->endC)) toend = 2; else toend = 4;
	ends = atend+toend;
	if (ends==2) { d = ats->endN|tos->endN; e = ats->endC|tos->endC; v = tos->endN-ats->endN; }
	if (ends==3) { d = ats->endC|tos->endN; e = ats->endN|tos->endC; v = tos->endN-ats->endC; }
	if (ends==4) { d = ats->endN|tos->endC; e = ats->endC|tos->endN; v = tos->endC-ats->endN; } 
	if (ends==5) { d = ats->endC|tos->endC; e = ats->endN|tos->endN; v = tos->endC-ats->endC; } 
	// Pr(d) NL
	if (d-e>bond-2 && d<bond+2) { // wrong ends are closer
		ats->spin();
		tos->spin();
	}
	if (d < bond) return;
	v.setVec(0.1);
	ats->move(v); tos->move(-v);
	// add link between basepairs
	at->parent->link[0].to = to->parent;
	at->parent->link[0].next = -1;
	FOR(i,ats->nbonds) {
		if (ats->bond[i].to ==0) continue;
		ats->bond[0].to = tos; // just for exempt()
		ats->bond[0].type = 1; // need to fix for draw()
		break;
	}
}

void setLinks ( float bondlen[][20] )
{ // preset link lengths for protein (bondlen[0]) and nucleic (bondlen[1])
float	bondCA = Data::bondCA,
	bondPP = Data::bondPP,
	distBP = 3.0;
	// Protein (bondlen[0]..)
	// SSE link lengths
	bondlen[0][0] = 1.7; bondlen[0][1] = 2.0;		// alpha SSE level
	bondlen[0][2] = 2.3; bondlen[0][3] = 2.6;		// alpha SSE level
	bondlen[0][4] = 2.9; bondlen[0][5] = 3.2;		// alpha SSE level
	bondlen[0][6] = bondlen[0][7] = bondCA*4.7/3.8;		// beta  SSE level
	bondlen[0][8] = bondlen[0][9] = bondCA*4.7/3.8;		// beta  SSE level
	// residue link length
	bondlen[0][10] = bondCA*5.1/3.8; 			// alpha CA (i->i+3)
	bondlen[0][11] = bondCA*6.2/3.8;			// alpha CA (i->i+4)
	bondlen[0][12] = bondlen[0][13] = bondCA*6.9/3.8;	// intra beta CA (i->i+/-2)
	bondlen[0][14] = bondlen[0][15] = bondCA*4.7/3.8;	// inter beta CA level
	bondlen[0][16] = Data::distH; // H+ on CA
	// RNA (bondlen[1]..)
	// base-pair link length
	distBP *= bondPP;
	bondlen[1][10] = distBP;			// RNA ladder
	bondlen[1][11] = distBP;			// 
}

void linker ( Cell *cell )
{
int	protein = 0;
float	bond, kick;
float	bondlen[4][20];
int	ctype,csort, ptype,psort;
int	nlinks = cell->nlinks,
	level = cell->level,
	model = cell->model;
Data	*param = Data::model+model;
float	basic = 1.5 * param->bumps[level];	// default link length is bump limit * 1.5
	if (basic<0) basic = -basic;
	moltype = param->moltype;
	subtype = param->subtype;
	if (moltype==0 && subtype==1) protein = 1;
	kick = LINK;
	if (cell==0) return;
	if (cell->empty) return;
//	if (level && cell->parent->solid > 0) return; 
	if (cell->link) { float d;
		ctype = cell->type; csort = cell->sort;
		ptype = cell->parent->type; psort = cell->parent->sort;
		setLinks(bondlen);
		FOR(i,nlinks) { Cell *link = cell->link[i].to; int type, brittle; float snap;
			if (link==0) continue;
			if (link->empty) continue;
		  	type = cell->link[i].type;
		  	snap = 0.01*(float)cell->link[i].next; // %stretch allowed before breaking
			if (snap < 0.0) brittle = 0; else brittle = 1; // -ve snap = can't break
			if (type && moltype<2) bond = bondlen[moltype][type]; else bond = basic;
			snap = bond*(snap+1.0);
			d = cell->xyz | link->xyz;
			if (brittle && d>snap) { // break an over stretched link
				link = 0;
				cell->link[i].to = 0;
				continue;
			}
			if (ptype==2 && psort==0) { // in loops, refine only bad links
				if (link==0) break;
				if (d > bond*1.5) {
					part2cells(cell, link, bond, kick);
				}
				if (d < bond*0.5) {
					part2cells(cell, link, bond, kick);
				}
				continue;
			}
			if (ctype==2 && csort) // in SSEs move along close approach not centres
			{ Vec	AonB = cell->xyz.vec_on_line(link->endN,link->endC),
				BonA = link->xyz.vec_on_line(cell->endN,cell->endC),
				midA = BonA & cell->xyz,
				midB = AonB & link->xyz,
				newA = midA, newB = midB;
				separate(newA,newB,bond,kick);
				cell->move(newA-midA);
				link->move(newB-midB);
				continue;
			}
			if (link==0) break;
			if (moltype==3) bond = 0.5*basic; // half atomic bump length for cells
			if (moltype==1) bond = 2.5; // basepair distance
			part2cells(cell, link, bond, kick);
		}
		// type 2 = SSE tube + sort 2 = beta (links 2 and 3 are between strands)
		if (protein && ptype==2 && psort==2) fixSheet(cell);
	}
	FOR(i,cell->kids) { //float far = -cell->far/(float)data[3];
		linker(cell->child[i]);
	}
}

void fixSheet ( Cell *b )
{
Cell	*a, *c;
Vec	s[3][3], oldcent, newcent;
Vec	A,B,C, x,y,z, mid, up, to, shift, told;
int	level = b->level, beta = 2;
Data	*param = Data::model+b->model;
float	size = param->bonds[level]+param->sizes[level];
float	bond = 3.8, span = 6.8, link = 4.7, drop = 1.5;
float	give = size*0.1, kick = 0.2, fix = 0.1, d;
int	i, j, flip0, flip2, inbeta;
float	w, wa, wb, wc, twist = 0.3;
	kick = fix = 0.1;
	span *= size/bond;
	link *= size/bond;
	drop *= size/bond;
	bond = size;
	// check all exist, are on the same level and most in beta
	if (!b->link) return;
	a = b->link[2].to;
	c = b->link[3].to;
	if (!a || !b || !c) return;
	// link->next holds the weight for the link (set -ve to prevent snapping in linker)
	wa = -(float)b->link[2].next;
	wc = -(float)b->link[3].next;
	wb = wa+wc;
	w = 1-exp(-wb*wb*0.1);
	if (randf() > w) return;
	wa = fix*wa;
	wb = fix*wb;
	wc = fix*wc;
	if (wb > 1.0) wb = 1.0;
	if (a->parent->sort != beta || b->parent->sort != beta || c->parent->sort != beta) return;
	if (a->sis->level != level) return; if (a->bro->level != level) return;
	if (b->sis->level != level) return; if (b->bro->level != level) return;
	if (c->sis->level != level) return; if (c->bro->level != level) return;
	inbeta = 0;
	if (a->sis->parent->sort == beta) inbeta++; if (a->bro->parent->sort == beta) inbeta++;
	if (b->sis->parent->sort == beta) inbeta++; if (b->bro->parent->sort == beta) inbeta++;
	if (c->sis->parent->sort == beta) inbeta++; if (c->bro->parent->sort == beta) inbeta++;
	if (inbeta < 5) return;
	// check atoms are far enough apart
	if (abs(a->atom - c->atom) < 6) return; 
	if (abs(a->atom - b->atom) < 4) return; 
	if (abs(b->atom - c->atom) < 4) return; 
	// check all lengths are reasonable before refining geometry
	d = vdif(a->xyz,a->sis->xyz)-bond;
	if (d*d>give) { part2cells(a,a->sis,bond,kick); return; }
	d = vdif(a->xyz,a->bro->xyz)-bond;
	if (d*d>give) { part2cells(a,a->bro,bond,kick); return; }
	d = vdif(b->xyz,b->sis->xyz)-bond;
	if (d*d>give) { part2cells(b,b->sis,bond,kick); return; }
	d = vdif(b->xyz,b->bro->xyz)-bond;
	if (d*d>give) { part2cells(b,b->bro,bond,kick); return; }
	d = vdif(c->xyz,c->sis->xyz)-bond;
	if (d*d>give) { part2cells(c,c->sis,bond,kick); return; }
	d = vdif(c->xyz,c->bro->xyz)-bond;
	if (d*d>give) { part2cells(c,c->bro,bond,kick); return; }
	d = vdif(a->sis->xyz,a->bro->xyz)-span;
	if (d*d>give) { part2cells(a->sis,a->bro,span,kick); return; }
	d = vdif(b->sis->xyz,b->bro->xyz)-span;
	if (d*d>give) { part2cells(b->sis,b->bro,span,kick); return; }
	d = vdif(c->sis->xyz,c->bro->xyz)-span;
	if (d*d>give) { part2cells(c->sis,c->bro,span,kick); return; }
	d = vdif(a->xyz,b->xyz)-link;
	if (d*d>give) { part2cells(a,b,link,kick); return; }
	d = vdif(c->xyz,b->xyz)-link;
	if (d*d>give) { part2cells(c,b,link,kick); return; }
	// fill working array (s) with parallel strands (s[strand][residue])
	s[0][0]=a->sis->xyz; s[0][1]=a->xyz; s[0][2]=a->bro->xyz;
	s[1][0]=b->sis->xyz; s[1][1]=b->xyz; s[1][2]=b->bro->xyz;
	s[2][0]=c->sis->xyz; s[2][1]=c->xyz; s[2][2]=c->bro->xyz;
	x = s[0][0]-s[0][2]; y = s[1][0]-s[1][2]; z = s[2][0]-s[2][2];
	flip0 = flip2 = 0;
	if (x*y < 0.0) { flip0 = 1;
		s[0][0]=a->bro->xyz; s[0][1]=a->xyz; s[0][2]=a->sis->xyz;
	}
	if (z*y < 0.0) { flip2 = 1;
		s[2][0]=c->bro->xyz; s[2][1]=c->xyz; s[2][2]=c->sis->xyz;
	}
	oldcent.zero();
	for (i=0; i<3; i++) for (j=0; j<3; j++) oldcent += s[i][j];
	oldcent /= 9.0;
	// set average strand in A-B-C (outer strands have half weight each)
	A = s[0][0] & s[2][0]; A = s[1][0]&A;
	B = s[0][1] & s[2][1]; B = s[1][1]&B;
	C = s[0][2] & s[2][2]; C = s[1][2]&C;
	separate(A,B,bond,kick);
	separate(C,B,bond,kick);
	separate(A,C,span,kick);
	up = B-(A&C);		// pleat direction 
	x = A-C;		// chain direction along average strand
	to = s[2][1]-s[0][1];	// mid strand to strand direction
	y = to^x;
	if (y*up < 0.0) y = -y; // keep Y pointing up
	z = x^y;
	if (z*to < 0.0) z = -z; // keep Z pointing as given
	// the xyz frame has x=strand(2-->0), y=up(pleat), z=across(0-->2)
	x.setVec(span*0.5);
	y.setVec(drop);
	z.setVec(link);
	// set A-B-C as B+/-x with pleat y
	A = B+x; A -= y;
	C = B-x; C -= y;
	// set new positions as average +/- z
	s[0][0]=A-z; s[0][1]=B-z; s[0][2]=C-z;
	s[1][0]=A;   s[1][1]=B;   s[1][2]=C;
	s[2][0]=A+z; s[2][1]=B+z; s[2][2]=C+z;
	// rotate corners to twist
	s[0][0].set_rot(s[2][1],s[0][1], twist);
	s[0][2].set_rot(s[2][1],s[0][1], twist);
	s[2][0].set_rot(s[2][1],s[0][1],-twist);
	s[2][2].set_rot(s[2][1],s[0][1],-twist);
	newcent.zero();
	for (i=0; i<3; i++) for (j=0; j<3; j++) newcent += s[i][j];
	newcent /= 9.0;
	told = oldcent-newcent;	// restore old centre
	for (i=0; i<3; i++) for (j=0; j<3; j++) s[i][j] += told;
	// update sheet coordinates (swap if flipped)
	if (flip0) { x=s[0][0]; s[0][0]=s[0][2]; s[0][2]=x; }
	if (flip2) { x=s[2][0]; s[2][0]=s[2][2]; s[2][2]=x; }
	//
	a->sis->xyz += (s[0][0] - a->sis->xyz)*wa*fix;
	a->xyz      += (s[0][1] - a->xyz     )*wa;
	a->bro->xyz += (s[0][2] - a->bro->xyz)*wa*fix;
	//
	b->sis->xyz += (s[1][0] - b->sis->xyz)*wb*fix;
	b->xyz      += (s[1][1] - b->xyz     )*wb;
	b->bro->xyz += (s[1][2] - b->bro->xyz)*wb*fix;
	//
	c->sis->xyz += (s[2][0] - c->sis->xyz)*wc*fix;
	c->xyz      += (s[2][1] - c->xyz     )*wc;
	c->bro->xyz += (s[2][2] - c->bro->xyz)*wc*fix;
}

Vec rdisp ( Vec a, Vec b, float r )
{ // add a random displacement (r) to b, keeping fixed length a-b
Vec     e, f;
float   d = (a|b);
	e = b;
	e.set_disp(r);
	f = (e-a).getVec(d);
	return Vec(a+f);
}

float scoreFit ( Vec a, Vec b, Vec c, Vec d, Vec e, Vec f, float leng ) 
{
float   ef,be,cf, dif, ta,td, score, wgap=100.0, woff=10.0, wang=1.0;
        be = b|e; cf = c|f; // keep close to old positions
        ta = PI - angle(a,e,f); // bias to towards a
        td = PI - angle(e,f,d); // convex connection
        ef = e|f;
        dif = ef-leng; dif *= dif;      // keep to ideal length
        score = dif*wgap+(be+cf)*woff+(ta+td)*wang;
//if (ef>0.0) { Pr(ef) Pr(dif) Pr(be) Pr(cf) Pr(ta) Pr(td) Pr(score) Pr(log(score)) NL } 
        return score;
}

#define TRY 10

void fixHinge ( Cell *at, Bonds *bond, float leng )
{
Vec	bs[TRY], cs[TRY];
Vec	a,b,c,d,e,f,g,h, u,v,w, x,y,z, beste, bestf;
int	i,j, m,n,in, link;
float	ab,bc,cd,ad, ra,rd, pa,pd, qa,qd, p,q;
float	ef, be, cf, dif, score, best, t,ta,td, max = 1.0;
Cell	*to = bond->to;
	link = bond->link;
	if (link==1) { a=at->endC; b=at->endN; c=to->endN; d=to->endC; } // NN
	if (link==2) { a=at->endC; b=at->endN; c=to->endC; d=to->endN; } // NC
	if (link==3) { a=at->endN; b=at->endC; c=to->endN; d=to->endC; } // CN
	if (link==4) { a=at->endN; b=at->endC; c=to->endC; d=to->endN; } // CN
	ab = a|b; bc = b|c; cd = c|d; ad = a|d;
	if (ad > ab+leng+cd) { float dist; // make colinear
		ta = angle(b,a,d);
		if (ta > NOISE) {
			if (ta>max) ta = max;
			x = (b-a).norm();
			y = (d-a).norm();
			z = x^y; g = a+z;
			spinCell(at,a,g,-ta);
		}
		td = angle(c,d,a);
		if (td>NOISE) {
			if (td>max) td = max;
			x = (c-d).norm();
			y = (a-d).norm();
			z = x^y; g = d+z;
			spinCell(to,d,g,-td);
		}
		dist = leng+(ab+cd)*0.5; // half cell separation
		part2cells(at,to,dist,0.5);
		return;
	}
	t =  leng-bc; if (t<0.0) t = -t;
	t = sqrt(t);
	if (t<0.01) return;
	ta = t*sqrt(ab)*1.0;
	td = t*sqrt(cd)*1.0;
	m = n = 1;
	bs[0] = b; cs[0] = c; // 1st try (bs[0],cs[0]) = current position
	best = scoreFit(a,b,c,d,b,c,leng);
	FOR(i,100) { float len2 = leng*2.0;;
		if (n==TRY && m==TRY) break;
		e = rdisp(a,bs[n-1],ta);
		f = rdisp(d,cs[m-1],td);
		score = scoreFit(a,b,c,d,e,f,leng);
		if (score > best) continue;
		in = 0;
		if (n<TRY && (e|d)>cd) { // e is outside D
			bs[n] = e;
			n++; in++;
		}
		if (m<TRY && vdif(f,a)>ab) { // f is outside A
			cs[m] = f;
			m++; in++;
		}
		if (in==2) best = score;
	}
	beste = b; bestf = c;
	FOR(i,n) { FOR(j,m) {
		score = scoreFit(a,b,c,d,bs[i],cs[j],leng);
		if (score<best) {
			best = score; beste = bs[i]; bestf = cs[j];
		}
	} }
	if ((b|beste) > 0.001) { // spin cell at
		ta = angle(b,a,beste);
		if (ta > NOISE) {
			if (ta>max) ta = max;
			x = (b-a).norm();
			y = (beste-a).norm();
			z = x^y; g = a+z;
			spinCell(at,a,g,-ta);
		}
	}
	if ((c|bestf) > 0.001) { // spin cell to
		td = angle(c,d,bestf);
		if (td>NOISE) {
			if (td>max) td = max;
			x = (c-d).norm();
			y = (bestf-d).norm();
			z = x^y; g = d+z;
			spinCell(to,d,g,-td);
		}
	}
} 
