#include "util.hpp"
#include "geom.hpp"
#include "cell.hpp"
#include "data.hpp"

#define BOND 0.5 // 0.3
#define LINK 0.5 // 0.6

void fixSheet ( Cell* );

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
		if (moltype==0 && end==0) {  // PROTein chain (not on an end)
			if (local[level]) {  // using preset lengths (set by -ve size in param.model
				sislen = cell->prox.x;	// use specific sis length (x)
				brolen = cell->prox.z;	// use specific bro length (z)
			}
			if (sislen >  999.0) sislen = size;	// use generic bond length
			if (brolen >  999.0) brolen = size;	// dito for bro
		}
		if (moltype==1) { // Nucleic acid chain
			if (level==depth-1) {	// don't refine basepaired stack bond distance to a loop
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
	DO(i,cell->kids) bonder(cell->child[i]);
}

void setLinks ( float bondlen[][20] )
{ // preset link lengths for protein (bondlen[0]) and nucleic (bondlen[1])
float	bondCA = Data::bondCA,
	bondPP = Data::bondPP,
	distBP = 3.0;
	// Protein
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
	// RNA
	// base-pair link length
	distBP *= bondPP;
	bondlen[1][10] = distBP;			// RNA ladder
	bondlen[1][11] = distBP;			// 
}

void linker ( Cell *cell )
{
float	bond, kick;
float	bondlen[4][20];
int	ctype,csort, ptype,psort;
int	nlinks = cell->nlinks,
	level = cell->level,
	model = cell->model;
Data	*param = Data::model+model;
float	basic = param->bumps[level];	// default link length is bump limit
	moltype = param->moltype;
	kick = LINK;
	if (cell==0) return;
	if (cell->empty) return;
//	if (level && cell->parent->solid > 0) return; 
	if (cell->link) { float d;
		ctype = cell->type; csort = cell->sort;
		ptype = cell->parent->type; psort = cell->parent->sort;
		setLinks(bondlen);
		DO(i,nlinks) { Cell *link = cell->link[i].to; int type; float snap;
			if (link==0) continue;
			if (link->empty) continue;
		  	type = cell->link[i].type;
		  	snap = 0.01*(float)cell->link[i].next; // %stretch allowed before breaking
			if (type && moltype<2) bond = bondlen[moltype][type]; else bond = basic;
			snap = bond*(snap+1.0);
			d = cell->xyz | link->xyz;
			if (d>snap) { // break a very long link
				link = 0;
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
			part2cells(cell, link, bond, kick);
		}
		// type 2 = SSE tube + sort 2 = beta (links 2 and 3 are between strands)
		if (moltype==0 && ptype==2 && psort==2) fixSheet(cell);
	}
	DO(i,cell->kids) { //float far = -cell->far/(float)data[3];
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
float	give = size*0.1, kick = 0.2, fix = 0.5, d;
int	i, j, flip0, flip2, inbeta;
float	wa, wb, wc, twist = 0.3;
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
	wa = fix*(float)b->link[2].next;
	wc = fix*(float)b->link[3].next;
	wb = fix*(wa+wc);
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
	a->sis->xyz += (s[0][0] - a->sis->xyz)*wa*fix;
	a->xyz      += (s[0][1] - a->xyz     )*wa;
	a->bro->xyz += (s[0][2] - a->bro->xyz)*wa*fix;
	b->sis->xyz += (s[1][0] - b->sis->xyz)*wb*fix;
	b->xyz      += (s[1][1] - b->xyz     )*wb;
	b->bro->xyz += (s[1][2] - b->bro->xyz)*wb*fix;
	c->sis->xyz += (s[2][0] - c->sis->xyz)*wc*fix;
	c->xyz      += (s[2][1] - c->xyz     )*wc;
	c->bro->xyz += (s[2][2] - c->bro->xyz)*wc*fix;
}

/*
#include "main/common.h"

char	num[11];

void linker ( int *datain, Cells *worldin, Types*** typesin )
{
int	i,j, m, n=0;
char	bondid[N+1];
float	*bondlen[2];
char	**bondstr[2];
float	distBP = 3.0;
//	bondxxx[0|1]: 0=prot, 1=RNA
	data = datain; world = worldin; types = typesin;
	if (data[0] > 0 && data[5] == 999 ) return; // freeze during presort
	for (i=0; i<2; i++) {
		bondlen[i] = (float*)alloca(sizeof(float)*22);
		bondstr[i] = (char**)alloca(sizeof(char*)*N);
		for (j=0; j<N; j++) bondstr[i][j] = (char*)alloca(sizeof(char*)*111);
	}
	depth = data[2];
	strcpy(num,"0123456789");
	//bondCA = 0.1*(float)(sizes[depth]+bonds[depth]);
	// Protein
	strcpy(bondstr[0][depth-1],"10-11-12-13-14-15-20-21-22-23-");  // SSE links
	// bondlen[0] match          0  1  2  3  4  5  6  7  8  9    to bondstr key
	strcpy(bondstr[0][depth],"100-101-200-201-202-203-");	// residue links
	// bondlen[0] match        10  11  12  13  14  15            to bondstr key
	// RNA
	strcpy(bondstr[1][depth],"100-101-200-201-202-203-");	// base pair links
	// bondlen[0] match        16  17  18  19  20  21            to bondstr key
	// Protein
	// SSE link lengths
	bondlen[0][0] = 1.7; bondlen[0][1] = 2.0;		// alpha SSE level
	bondlen[0][2] = 2.3; bondlen[0][3] = 2.6;		// alpha SSE level
	bondlen[0][4] = 2.9; bondlen[0][5] = 3.2;		// alpha SSE level
	bondlen[0][6] = bondlen[0][7] = -bondCA*4.7/3.8;	// beta  SSE level (-ve = only bumps)
	bondlen[0][8] = bondlen[0][9] = -bondCA*4.7/3.8;	// beta  SSE level (-ve = only bumps)
	// residue link length
	bondlen[0][10] = bondCA*5.1/3.8; 			// alpha CA (i->i+3)
	bondlen[0][11] = bondCA*6.2/3.8;			// alpha CA (i->i+4)
	bondlen[0][12] = bondlen[0][13] = bondCA*6.9/3.8;	// intra beta CA (i->i+/-2)
	bondlen[0][14] = bondlen[0][15] = bondCA*4.7/3.8;	// inter beta CA level
	// RNA
	// base-pair link length
	distBP *= bondPP;
	bondlen[1][10] = distBP;			// RNA ladder
	bondlen[1][11] = distBP;			// 
	linkCell(world,0,bondid,bondlen,bondstr);
	bondCell(world,0,bondid,bondlen,bondstr);
}

linkCell ( Cells *cell, int level, char *bondid, float **bondlenin, char ***bondstrin )
{
float	d;
Vec	shift;
int	i, j, n, m, id;
Cells	*a, *b;
float	size, bond, basic, kick, len;
float	*bondlen;
char	**bondstr;
int	nbonds;
	if (cell==0) return;
	if (cell->empty) return;
//	if (level && cell->parent->solid > 0) return; 
	id = cell->id;
	n = cell->kids;
	m = cell->model*M*N+N;
	moltype = data[m];
	bondlen = bondlenin[moltype];
	bondstr = bondstrin[moltype];
	align = data+m+N*1; class = data+m+N*2;
        sizes = data+m+N*3; bumps = data+m+N*4; links = data+m+N*5; chain = data+m+N*6;
	kicks = data+m+N*7; keeps = data+m+N*8; repel = data+m+N*9; bonds = data+m+N*10;
	split = data+m+N*11; local = data+m+N*12; // flags local or generic bond length
	kick = 0.2; // WAS 0.1;
	size = 0.1*(float)(sizes[level]+bonds[level]);
	bondid[level] = num[cell->sort];
	basic = 0.1*(float)bumps[level]; // default link replaced later from types (if set)
	if (cell->link) { float d, error = 0.0;
		bondid[level+2] = '\0';
		for (i=0; i<links[level]; i++) { char *bids; int bid;
			if (cell->link[i]==0) continue;
			if (cell->link[i]->empty) continue;
			d = vdif(cell->xyz,cell->link[i]->xyz);
			if (d>30.0) { // break a very long link
				cell->link[i] = 0;
				continue;
			}
			bond = basic;
			bondid[level+1] = num[i];
			bondid[level+2] = '-';
			bondid[level+3] = '\0';
			bids = strstr(bondstr[level],bondid+depth-1);
		       	if (bids) {
				bid = (int)(bids-bondstr[level]);
				if (level==depth) bid = bid/4 + 10; else bid = bid/3;
				bond = bondlen[bid];
				if (bond < 0.0) {	// just check bumps
					if (vdif(cell->xyz,cell->link[i]->xyz) > bond) continue;
					bond = -bond;
				}
			} 
			if (cell->type==2 && cell->sort==0) { // don't refine loops much
				if (cell->link[i]==0) break;
				d = vdif(cell->xyz,cell->link[i]->xyz);
				if (d > bond*1.5) {
					part2cells(cell, cell->link[i], bond, kick, 1);
				}
				if (d < bond*0.5) {
					part2cells(cell, cell->link[i], bond, kick, 1);
				}
				continue;
			}
			if (cell->link[i]==0) break;
			//  stops reactions? bond += (float)(depth-level);
			part2cells(cell, cell->link[i], bond, kick, 0); // don't use knockon = 1);
			if (cell->type == 2)	// keep space between SSE ends (-kick = repel only)
			{ float dnn = vdif(cell->endN,cell->link[i]->endN),
				dnc = vdif(cell->endN,cell->link[i]->endC),
				dcn = vdif(cell->endC,cell->link[i]->endN),
				dcc = vdif(cell->endC,cell->link[i]->endC);
				if (dnn+dcc < dnc+dcn) {	// parallel 
					separate(&(cell->endN),&(cell->link[i]->endN),bond,-kick*0.5);
					separate(&(cell->endC),&(cell->link[i]->endC),bond,-kick*0.5);
				} else {			// antiparr
					separate(&(cell->endN),&(cell->link[i]->endC),bond,-kick*0.5);
					separate(&(cell->endC),&(cell->link[i]->endN),bond,-kick*0.5);
				}
			}
		}
		bondid[level+1] = '\0';
		// type 2 = SSE tube + sort 2 = beta and links 2 and 3 are between strands
		if (level && cell->parent->type==2 && cell->parent->sort==2) {
			fixSheet(cell->link[2], cell, cell->link[3]);
		}
	}
	for (i=0; i<n; i++) { //float far = -cell->far/(float)data[3];
		linkCell(cell->child[i],level+1,bondid,bondlenin,bondstrin);
	}
}

bondCell ( Cells *cell, int level, char *bondid, float **bondlenin, char ***bondstrin )
{
float	d;
Vec	shift;
int	i, j, n, m, id;
Cells	*a, *b;
float	size, bond, basic, kick, len;
float	*bondlen;
char	**bondstr;
int	nbonds;
	if (cell==0) return;
	if (cell->empty) return;
	id = cell->id;
	n = cell->kids;
//	if (level && cell->parent->solid > 0) return; 
	m = cell->model*M*N+N;
	moltype = data[m];
	align = data+m+N*1; class = data+m+N*2;
        sizes = data+m+N*3; bumps = data+m+N*4; links = data+m+N*5; chain = data+m+N*6;
	kicks = data+m+N*7; keeps = data+m+N*8; repel = data+m+N*9; bonds = data+m+N*10;
	split = data+m+N*11; local = data+m+N*12; // flags local or generic bond length
	len = size = 0.1*(float)(sizes[level]+bonds[level]);
	kick = 0.3;
	DO(j,1) {
		if (moltype < 2) {
                	if (cell->bond && cell->bond[0].type ) { // use HINGE data for protein or nucleic
                        	for (i=0; i<5; i++) { Bonds *hinge = cell->bond+i;
                                	if (hinge->to == 0) break;
                                	fixHinge(cell,hinge);
                        	}
				break;
			}
                }
		if (moltype==0) { // PROTein chain
			if (local[level]) { // using preset lengths (set by -ve size in param.model
				len = cell->prox.x;	// use specific sis length (x)
			}
			if (chain[level]==0 || (cell->sis->level!=cell->level)) break;
			if (len > 999.0) len = size;  // use generic bond length
			part2cells(cell,cell->sis,len,kick,0);
			break;
		}
		if (moltype==1) { // nucleic
			if (chain[level]==0 || (cell->sis->level!=cell->level)) break;
			if (level==depth) { // atomic
				part2cells(cell,cell->sis,size,kick,0);
			}
			if (level==depth-1 && (int)cell->endC.z!=1234 && (int)cell->endN.z!=1234   // z = 1234
			         && (int)cell->sis->endC.z!=1234 && (int)cell->sis->endN.z!=1234 ) // for unset
			{ char	stemcel, stemsis; // SSE, endC is near termini
			  float gap, disp; Vec move;
				if (cell->sort==1) stemcel = 1; else stemcel = 0;
				if (cell->sis->sort==1) stemsis = 1; else stemsis = 0;
				if ( stemcel &&  stemsis) vsub(cell->endC,cell->sis->endC, &move);
				if ( stemcel && !stemsis) vsub(cell->endC,cell->sis->endC, &move);
				if (!stemcel &&  stemsis) vsub(cell->endN,cell->sis->endC, &move);
				if (!stemcel && !stemsis) vsub(cell->endC,cell->sis->endC, &move);
				gap = vmod(move); disp = size-gap;
				vnorm(&move); vmul(&move,disp*kick*0.5);
				moveCell(cell->sis,move,-1);
				moveCell(cell,move, 1);
			}
			if (level<depth-1) { // other (as polymer)
				part2cells(cell,cell->sis,size,kick,0);
			}
			break;
		}
		if (moltype==2) { // CHEMical bonds
			nbonds = cell->nbonds;
			for (i=0; i<nbonds; i++) {
				if (cell->bond[i].to == 0) continue; // cyclic bonds at end
				part2cells(cell,cell->bond[i].to,size,kick,0);
			}
			break;
		}
	}
	for (i=0; i<n; i++) { // float far = -cell->far/(float)data[3];
		bondCell(cell->child[i],level+1,bondid,bondlenin,bondstrin);
	}
}

rdisp ( Vec a, Vec b, Vec *c, float r )
{ // add a random displacement (r) to b, keeping fixed length a-b
Vec	e, f;
float	d = vdif(a,b);
	vcopy(b,&e);
	vradd(&e,r);
	vsub(e,a,&f);
	vnorm(&f);
	vmul(&f,d);
	vadd(a,f,c);
}

float scoreFit ( Vec a, Vec b, Vec c, Vec d, Vec e, Vec f, float leng) 
{
float	ef,be,cf, dif, ta,td, score, wgap=100.0, woff=10.0, wang=1.0;
	be = vdif(b,e); cf = vdif(c,f); // keep close to old positions
	ta = PI - angle(a,e,f); // bias to towards a
	td = PI - angle(e,f,d); // convex connection
	ef = vdif(e,f);
	dif = ef-leng; dif *= dif;	// keep to ideal length
	score = dif*wgap+(be+cf)*woff+(ta+td)*wang;
//if (ef>0.0) { Pr(ef) Pr(dif) Pr(be) Pr(cf) Pr(ta) Pr(td) Pr(score) Pr(log(score)) NL } 
	return score;
}

#define TRY 10

fixHinge (Cells *at, Bonds *bond)
{
Vec	bs[TRY], cs[TRY];
Vec	a,b,c,d,e,f,g,h, u,v,w, x,y,z, beste, bestf;
int	i,j, m,n,in, link;
float	leng, ab,bc,cd,ad, ra,rd, pa,pd, qa,qd, p,q;
float	ef, be, cf, dif, score, best, t,ta,td, max = 1.0;
Cells	*to;
	to = bond->to;
	link = bond->next;
	leng = 0.1*(float)bond->type; 
	if (link==1) { vcopy(at->endC,&a); vcopy(at->endN,&b); vcopy(to->endN,&c); vcopy(to->endC,&d); } // NN
	if (link==2) { vcopy(at->endC,&a); vcopy(at->endN,&b); vcopy(to->endC,&c); vcopy(to->endN,&d); } // NC
	if (link==3) { vcopy(at->endN,&a); vcopy(at->endC,&b); vcopy(to->endN,&c); vcopy(to->endC,&d); } // CN
	if (link==4) { vcopy(at->endN,&a); vcopy(at->endC,&b); vcopy(to->endC,&c); vcopy(to->endN,&d); } // CN
	ab = vdif(a,b); bc = vdif(b,c); cd = vdif(c,d); ad = vdif(a,d);
	if (ad > ab+leng+cd) { float dist; // make colinear
		ta = angle(b,a,d);
		if (ta > NOISE) {
			if (ta>max) ta = max;
			vsub(b,a,&x); vnorm(&x);
			vsub(d,a,&y); vnorm(&y);
			vprod(x,y,&z); vadd(a,z,&g);
			spinCell(at,a,g,-ta);
		}
		td = angle(c,d,a);
		if (td>NOISE) {
			if (td>max) td = max;
			vsub(c,d,&x); vnorm(&x);
			vsub(a,d,&y); vnorm(&y);
			vprod(x,y,&z); vadd(d,z,&g);
			spinCell(to,d,g,-td);
		}
		dist = leng+(ab+cd)*0.5; // half cell separation
		part2cells(at,to,dist,0.5,0);
		return;
	}
	t =  leng-bc; if (t<0.0) t = -t;
	t = sqrt(t);
	if (t<0.01) return;
	ta = t*sqrt(ab)*1.0;
	td = t*sqrt(cd)*1.0;
	m = n = 1;
	vcopy(b,bs); vcopy(c,cs); // 1st try (bs[0],cs[0]) = current position
	best = scoreFit(a,b,c,d,b,c,leng);
	DO(i,1000) { float len2 = leng*2.0;;
		if (n==TRY && m==TRY) break;
		rdisp(a,bs[n-1],&e,ta);
		rdisp(d,cs[m-1],&f,td);
		score = scoreFit(a,b,c,d,e,f,leng);
		if (score > best) continue;
		in = 0;
		if (n<TRY && vdif(e,d)>cd) { // e is outside D
			vcopy(e,bs+n);
			n++; in++;
		}
		if (m<TRY && vdif(f,a)>ab) { // f is outside A
			vcopy(f,cs+m);
			m++; in++;
		}
		if (in==2) best = score;
	}
	vcopy(b,&beste); vcopy(c,&bestf);
	DO(i,n) { DO(j,m) {
		score = scoreFit(a,b,c,d,bs[i],cs[j],leng);
		if (score<best) {
			best = score; vcopy(bs[i],&beste); vcopy(cs[j],&bestf);
		}
	}	}
	if (vdif(b,beste) > 0.001) { // spin cell at
		ta = angle(b,a,beste);
		if (ta > NOISE) {
			if (ta>max) ta = max;
			vsub(b,a,&x); vnorm(&x);
			vsub(beste,a,&y); vnorm(&y);
			vprod(x,y,&z); vadd(a,z,&g);
			spinCell(at,a,g,-ta);
		}
	}
	if (vdif(c,bestf) > 0.001) { // spin cell to
		td = angle(c,d,bestf);
		if (td>NOISE) {
			if (td>max) td = max;
			vsub(c,d,&x); vnorm(&x);
			vsub(bestf,d,&y); vnorm(&y);
			vprod(x,y,&z); vadd(d,z,&g);
			spinCell(to,d,g,-td);
		}
	}
} 

fixSheet (Cells *a, Cells *b, Cells *c)
{
Vec	s[3][3], oldcent, newcent;
Vec	A,B,C, x,y,z, mid, up, to, shift, told;
float	size = 0.1*(float)(bonds[b->level]+sizes[b->level]);
float	bond = 3.8, span = 6.8, link = 4.7, drop = 1.5;
float	give = size*0.1, kick = 0.2, fix = 0.5, d;
int	i, j, flip0, flip2, inbeta;
float	twist = 0.3;
int	level = b->level, beta = 2;
	span = size*span/bond;
	link = size*link/bond;
	drop = size*drop/bond;
	bond = size;
	kick = 0.1; fix = 0.01;
	// check all exist, are on the same level and most in beta
	if (!a || !b || !c) return;
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
	if (d*d>give) { separate(&(a->xyz),&(a->sis->xyz),bond,kick); return; }
	d = vdif(a->xyz,a->bro->xyz)-bond;
	if (d*d>give) { separate(&(a->xyz),&(a->bro->xyz),bond,kick); return; }
	d = vdif(b->xyz,b->sis->xyz)-bond;
	if (d*d>give) { separate(&(b->xyz),&(b->sis->xyz),bond,kick); return; }
	d = vdif(b->xyz,b->bro->xyz)-bond;
	if (d*d>give) { separate(&(b->xyz),&(b->bro->xyz),bond,kick); return; }
	d = vdif(c->xyz,c->sis->xyz)-bond;
	if (d*d>give) { separate(&(c->xyz),&(c->sis->xyz),bond,kick); return; }
	d = vdif(c->xyz,c->bro->xyz)-bond;
	if (d*d>give) { separate(&(c->xyz),&(c->bro->xyz),bond,kick); return; }
	d = vdif(a->sis->xyz,a->bro->xyz)-span;
	if (d*d>give) { separate(&(a->sis->xyz),&(a->bro->xyz),span,kick); return; }
	d = vdif(b->sis->xyz,b->bro->xyz)-span;
	if (d*d>give) { separate(&(b->sis->xyz),&(b->bro->xyz),span,kick); return; }
	d = vdif(c->sis->xyz,c->bro->xyz)-span;
	if (d*d>give) { separate(&(c->sis->xyz),&(c->bro->xyz),span,kick); return; }
	d = vdif(a->xyz,b->xyz)-link;
	if (d*d>give) { separate(&(a->xyz),&(b->xyz),link,kick); return; }
	d = vdif(c->xyz,b->xyz)-link;
	if (d*d>give) { separate(&(c->xyz),&(b->xyz),link,kick); return; }
	// fill working array (s) with parallel strands
	vcopy(a->sis->xyz,&(s[0][0])); vcopy(a->xyz,&(s[0][1])); vcopy(a->bro->xyz, &(s[0][2]));
	vcopy(b->sis->xyz,&(s[1][0])); vcopy(b->xyz,&(s[1][1])); vcopy(b->bro->xyz, &(s[1][2]));
	vcopy(c->sis->xyz,&(s[2][0])); vcopy(c->xyz,&(s[2][1])); vcopy(c->bro->xyz, &(s[2][2]));
	vsub(s[0][0],s[0][2], &x); vsub(s[1][0],s[1][2], &y); vsub(s[2][0],s[2][2], &z);
	flip0 = flip2 = 0;
	if (vdot(x,y) < 0.0) {
		vcopy(a->bro->xyz,&(s[0][0])); vcopy(a->xyz,&(s[0][1])); vcopy(a->sis->xyz, &(s[0][2]));
		flip0 = 1;
	}
	if (vdot(z,y) < 0.0) {
		vcopy(c->bro->xyz,&(s[2][0])); vcopy(c->xyz,&(s[2][1])); vcopy(c->sis->xyz, &(s[2][2]));
		flip2 = 1;
	}
	vinit(&oldcent);
	for (i=0; i<3; i++) for (j=0; j<3; j++) vsum(s[i][j], &oldcent);
	vdiv(&oldcent, 9.0);
	// set average strand in A-B-C (outer strands half weight)
	vave(s[0][0], s[2][0], &A); vave(s[1][0], A, &A);
	vave(s[0][1], s[2][1], &B); vave(s[1][1], B, &B);
	vave(s[0][2], s[2][2], &C); vave(s[1][2], C, &C);
	separate(&A,&B,bond,kick);
	separate(&C,&B,bond,kick);
	separate(&A,&C,span,kick);
	vave(A,C, &mid);
	vsub(B,mid, &up);		// pleat direction 
	vsub(A, C, &x);			// chain direction
	vsub(s[2][1], s[0][1], &to);	// strand to strand direction
	vprod(to,x, &y);
	if (vdot(y,up) < 0.0) vmul(&y,-1.0); // keep Y pointing up
	vprod(x,y, &z);
	if (vdot(z,to) < 0.0) vmul(&z,-1.0); // keep Z pointing as given
	// the xyz frame has x=strand, y=up(pleat), Z=2-->0
	vnorm(&x); vmul(&x,span*0.5);
	vnorm(&y); vmul(&y,drop);
	vnorm(&z); vmul(&z,link);	// set new positions as average +/- z
	vadd(B,x,&A); vsub(A,y,&A);
	vsub(B,x,&C); vsub(C,y,&C);
	vsub(A,z,&(s[0][0])); vsub(B,z,&(s[0][1])); vsub(C,z,&(s[0][2]));
	vcopy(A, &(s[1][0])); vcopy(B, &(s[1][1])); vcopy(C, &(s[1][2]));
	vadd(A,z,&(s[2][0])); vadd(B,z,&(s[2][1])); vadd(C,z,&(s[2][2]));
	// rotate corners to twist
	rotate(s[2][1],s[0][1],&(s[0][0]), twist);
	rotate(s[2][1],s[0][1],&(s[0][2]), twist);
	rotate(s[2][1],s[0][1],&(s[2][0]),-twist);
	rotate(s[2][1],s[0][1],&(s[2][2]),-twist);
	vinit(&newcent);
	for (i=0; i<3; i++) for (j=0; j<3; j++) vsum(s[i][j], &newcent);
	vdiv(&newcent, 9.0);
	vsub(oldcent,newcent, &told);	// restore old centre
	for (i=0; i<3; i++) for (j=0; j<3; j++) vsum(told, &(s[i][j]));
	// update sheet coordinates (swap if flipped)
	if (flip0) { vcopy(s[0][0],&x); vcopy(s[0][2],&(s[0][0])), vcopy(x,&(s[0][2])); }
	if (flip2) { vcopy(s[2][0],&x); vcopy(s[2][2],&(s[2][0])), vcopy(x,&(s[2][2])); }
	vsub(s[0][0], a->sis->xyz, &shift); vmul(&shift,fix); vsum(shift, &(a->sis->xyz));
	vsub(s[0][1], a->xyz,      &shift); vmul(&shift,fix); vsum(shift, &(a->xyz));
	vsub(s[0][2], a->bro->xyz, &shift); vmul(&shift,fix); vsum(shift, &(a->bro->xyz));
	vsub(s[1][0], b->sis->xyz, &shift); vmul(&shift,fix); vsum(shift, &(b->sis->xyz));
	vsub(s[1][1], b->xyz,      &shift); vmul(&shift,fix); vsum(shift, &(b->xyz));
	vsub(s[1][2], b->bro->xyz, &shift); vmul(&shift,fix); vsum(shift, &(b->bro->xyz));
	vsub(s[2][0], c->sis->xyz, &shift); vmul(&shift,fix); vsum(shift, &(c->sis->xyz));
	vsub(s[2][1], c->xyz,      &shift); vmul(&shift,fix); vsum(shift, &(c->xyz));
	vsub(s[2][2], c->bro->xyz, &shift); vmul(&shift,fix); vsum(shift, &(c->bro->xyz));
}
*/
