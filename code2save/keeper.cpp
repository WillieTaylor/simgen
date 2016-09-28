#include "util.hpp"
#include "geom.hpp"
#include "cell.hpp"
#include "data.hpp"

int  packBall ( Cell*, Cell*, float );
int  packTube ( Cell*, Cell*, float );
void packBase ( Cell*, Cell*, float );
int  unpackTube ( Cell*, Cell*, float );
float inEgg ( Vec, Vec, Vec, float, Vec& );
float inEgg ( Vec, Vec, Vec, float );
float inEgg ( Vec, Seg, float );
int packEgg ( Cell*, Cell*, float );

float margin = 0.9, margout = 1.0/margin;

void center ()
{
Cell	*world = Cell::world;
	FOR(i,world->kids) { Cell *child = world->child[i];
                if (child->empty) continue;
		child->group();
	}
}

void keeper ( Cell *cell )
{
int	m = cell->model, n = cell->kids,
	moltype = Data::model[m].moltype,
	subtype = Data::model[m].subtype,
	dna = moltype*subtype, run = 1;
	if ( randf() < 0.5 ) run = -1;
	FOR(i,n) // push the children inside
	{ Cell *child;
	  int	type, sort, lost = 0;
	  float s, strict, balance;
		if ( run > 0 ) {	// forward order
			child = cell->child[i];
		} else {		// backwards
			child = cell->child[n-i-1];
		}
                if (child->empty) continue;
	  	type = cell->type; sort = cell->sort;
/* ?
		if (chain[level] > 0 && chain[level+1] > 0) {   // in a chain of chains (don't push ends in)
			if (child == cell->starts || child == cell->finish) continue;
		}
*/
		strict = 0.1;  /////////// WAS 1.0
		switch (type) {
			case 0 :	// dummy sphere
				packBall(cell,child,1.0);
			break;
			case 1 :	// simple sphere
				packBall(cell,child,strict);
			break;
			case 2 :	// tube
				if (cell->level==depth-1 && sort==0) s = 0.1*strict; else s = strict;
				if (moltype==1) { // nucleic
					packBase(cell,child,s);
				} else {	  // protein
					packTube(cell,child,s);
				}
			break;
			case 3 :	// Ellipsoid
				if (cell->sort == E/2) {	// same as a sphere
					packBall(cell,child,1.0);
				} else {
					packEgg(cell,child,1.0);
				}
			break;
		}
	}
	FOR(i,cell->kids) keeper(cell->child[i]);
}

void packBase ( Cell *cell, Cell* child, float strict ) {
// <cell> = 'basepair secondary structure', child = base (P), <pa> = segment (domain)
Cell	*pa = cell->parent;
Vec	mid, shift;
	if (pa->type != 2 ) return;
	if (child->level != depth) return;
	// confine atoms to level-2 tube (DNA segment is at domain level)
	if (cell->sort == 0) { // don't confine loops (expel)
		unpackTube(pa,child,strict); // NB: children packed into parent level tube
	} else {
		packTube(pa,child,strict); // NB: children packed into parent level tube
	}
	if (cell->sort == 0) {     // realign loop tube with ends of the chain
		shift = (cell->junior->xyz - cell->senior->xyz).norm();
		cell->endN = cell->xyz - shift;
		cell->endC = cell->xyz + shift;
		return;
	}
	// reset basepair tube ends
	part2cells(cell->child[0],cell->child[1],2.7,0.1);
	cell->xyz = cell->child[0]->xyz & cell->child[1]->xyz;
	shift = (cell->child[1]->xyz - cell->child[0]->xyz)*0.7/2.0;
	cell->endN = cell->xyz - shift;
	cell->endC = cell->xyz + shift;
	// shift basepair tube centre towards axis
	mid = (cell->xyz).vec_on_line(pa->endN,pa->endC);
	shift = mid - cell->xyz;
	cell->move(shift*0.1);
}

int packBall ( Cell* cell, Cell* child, float strict )
{
int	shell;
Vec	shift;
Data	*model = Data::model+cell->model;
float	radius = model->sizes[cell->level]*0.5,
	push   = model->keeps[cell->level]*strict, d;
	if (cell->level==Data::depth) return 0;
	if ((int)cell->endN.z == 1234) return 0;	// ends not set yet
	if ((int)cell->endC.z == 1234) return 0;	// ends not set yet
	if (model->moltype < 3) radius -= model->sizes[child->level]*0.5; // keep totally inside
	if (push < -NOISE) shell = 1;			// confine to shell (+/-margin)
		 else	   shell = 0;			// confine inside
	shift = cell->xyz - child->xyz;			// shift from child to zero
	d = shift.len();
	if (shell) {
		if (d>radius*margin && d<radius*margout) return 0;
		if (d > radius) push = -push;		// outside the sphere
	} else {
		if (d < radius*margout) return 0;	// within the sphere
	}
	shift *= push;
	if (shift.sqr()<NOISE) return 0;
	child->move(shift);
	return 1;
}

int unpackTube ( Cell *cell, Cell *child, float strict )
{
int	in, shell = 1, closed = 1; //Data::tubecaps;
Vec	shift, over, temp;
float	dlast, dnext, dif, d, pale;
Data	*model = Data::model+cell->model;
float	radius = model->sizes[cell->level]*0.5,
	push   = model->keeps[cell->level]*strict;
int	moltype = model->moltype,
	subtype = model->subtype,
	protein = 0;
	if (moltype==0 && subtype==1) protein = 1;
	if (cell->level==Data::depth) return 0;
	if ((int)cell->endN.z == 1234) return 0; // ends not set yet
	if ((int)cell->endC.z == 1234) return 0; // ends not set yet
	if (strict < 0) strict = -strict;
	pale = radius*1.5;
	d = cell->xyz|child->xyz;
	if (d > pale*2.0) return 0;
	if (d < pale) {
		push = (pale-d)*strict;
		shift = (child->xyz - cell->xyz).norm()*push;
		child->move(shift);
	}
	pale = radius*1.1;
	d = cell->endN|child->xyz;
	if (d < pale) {
		push = (pale-d)*strict;
		shift = (child->endN - cell->xyz).norm()*push;
		child->move(shift);
	}
	d = cell->endC|child->xyz;
	if (d < pale) {
		push = (pale-d)*strict;
		shift = (child->endC - cell->xyz).norm()*push;
		child->move(shift);
	}
}

int packTube ( Cell *cell, Cell *child, float strict )
{
int	in, shell = 1, closed = 1; //Data::tubecaps;
Vec	shift, over, temp;
float	dlast, dnext, dif, d, pale;
Data	*model = Data::model+cell->model;
float	radius = model->sizes[cell->level]*0.5,
	push   = model->keeps[cell->level]*strict;
int	moltype = model->moltype,
	subtype = model->subtype,
	protein = 0;
	if (moltype==0 && subtype==1) protein = 1;
	if (cell->level==Data::depth) return 0;
	if ((int)cell->endN.z == 1234) return 0; // ends not set yet
	if ((int)cell->endC.z == 1234) return 0; // ends not set yet
	if (protein && cell->level==depth-1) { // apply SSE factor
		if (cell->sort==0) radius *= loopTHIC;
		if (cell->sort==1) radius *= alphTHIC;
		if (cell->sort==2) radius *= betaTHIC;
	}
//	radius -= model->sizes[child->level]*0.5; // keep totally inside
	pale = radius;
	if (push < -NOISE) {	// -ve push = flag to confine inside
		shell = 0; push = -push;
	}
	if (protein && cell->sort==0) {
		shell = closed = 0;	// protein loops always have open ends and no shell
	}
	if (vdif(cell->endN,cell->endC) < NOISE) return 0;
	if (closed) {
		dlast = cell->endN|child->xyz;
		dnext = cell->endC|child->xyz;
		if (!shell) { // make quick pre-check at ends
			if (dlast < pale) return 0;
			if (dnext < pale) return 0;
		}
	}
	in   = child->xyz.vec_in_seg (cell->endN, cell->endC); // in the line segment
	d    = child->xyz.vec_to_line(cell->endN, cell->endC); // dist to extended line
	over = child->xyz.vec_on_line(cell->endN, cell->endC); // image on extended line
	shift = over - child->xyz;
	dif = fabs(pale-d);
	push *= dif;
	if (closed) { // confine within segment
		if (in) { // over line segment
			if (d < pale*margout) { // inside the tube
				if (shell==0) return 0;	// inside is OK
				if (d < pale*margin) {	// too deep inside
					push = -push;	// push out
				} else { return 0; }	// in the shell margin
			}
		} else { // over an end-cap (so find which one)
			if (dlast < dnext) {		// over endN
				shift = cell->endN - child->xyz;
				d = cell->endN | child->xyz;
			} else{				// over endC
				shift = cell->endC - child->xyz;
				d = cell->endC | child->xyz;
			}
			dif = fabs(pale-d);
			push *= dif;
			if (d < pale*margout) { // inside the cap
				if (shell==0) return 0;	// inside is OK
				if ( d < pale*margin) { // too deep inside
					dif = fabs(pale-d);
					push = -push;	// push out
				} else { return 0; }	// in shell margin
			}
		}
	} else {	// confine within extended tube
		if (d < pale*margout) { // inside the tube
			if (shell==0) return 0;	// inside is OK
			if (d < pale*margin) {	// too deep inside
				push = -push;	// push out
			} else { return 0; }
		}
	}
	shift.setVec(push);
	if (shift.sqr()<NOISE) return 0;
	child->move(shift);
	return 1;
}

float inEgg ( Vec cell, Vec endN, Vec endC, float len, Vec &move )
// returns a value that is +ve outside and -ve inside the ellipsoid (not surface distance)
// cell = point, endN,endC = poles, len = diameter across minor axes, move = returned shift
{
float	dB, rB, dA, rA;
Vec	p,q,r, cent, axis;
float	a, b, c, d, e, f, g;
float	v,w,x,y,z;
int	apply = 1;
	if (len<0) { len = -len; apply = 0; }
	cent = endN & endC;	// set centre
	axis = endC - endN;	// axis between poles
	dA = endN|endC;		// dist between poles
	dB = len;		// diameter across minor axis
	rA = dA*0.5;		// major (rA) and...
	rB = dB*0.5;		// minor (rB) semi-axis lengths
	x = y = 1.0; z = dA/len;
	// bubble sort axes to x>y>z
	if (x<y) { w=x; x=y; y=w; }
	if (y<z) { w=y; y=z; z=w; }
	if (x<y) { w=x; x=y; y=w; }
	v = y/x; w = z/x;
	d = rA;	// major semi-axis length
	a = 1.0-v; b = 1.0-w; // fractional axis lengths relative to X
	if (a*a+b*b < 0.04) { // sphere-ish (both minor axes over 0.8-ish)
		g = 3.0*d/(1.0+v+w);	// mean sphere radius
		c = cell|cent;
		if (apply) move = (cent-cell).norm();	// move towards centre
		return c-g;
	}
	if (v > 0.5+0.5*w) {	// oblate
		g = d;
		d = d/w*(1.0+v)*0.5; // force radial symmetry about Z (min)
		p = cell-cent; q = axis; r = p^q;
		axis = r^q;	// axis now lies on cell-cent-pole plane
	} else {		//prolate
		g = d*0.5*(v+w);	// force radial symmetry about X (max)
	}
/*
                             * cell (external point)
                     ..-+-../
                .    d /|  /    .           d = focus pole distance
              .       /b| /       .         b = minor axis length/2
             :       /  |/c        :        c = foucs cent distance
        endN |------x---+---x------| endC   x = focii on major axis
                       cent                 a = major axis length/2
	length of focus1-surface-focus2 path:
	at major axis = c+a+(a-c) = 2a
	at minor axis = 2d, dd = bb+cc
	since the paths are equal: d = a;
	so	cc = dd-bb = aa-bb
	and 	c = sqrt(aa-bb)
*/
	c = sqrt(d*d-g*g);	// distance from centre to focii
	axis.setVec(c);		// set axis to length c
	p = cent+axis;		// up axis to focus1
	e = cell|p;		// dist to focus1
	q = cent-axis;		// down axis to focus2
	f = cell|q;		// dist of focus2
	if (apply) { // move is p*f+q*e (swap e/f weights to bias towards nearest focus)
		p -= cell;	// vect to focus1
		q -= cell;	// vect to focus2
		p.setVec(f);
		q.setVec(e);
		move = (p+q).norm();
	}
	return 1.4*(0.5*(e+f)-d);
}

float inEgg ( Vec cell, Vec endN, Vec endC, float len )
{ // just return the distance (-ve len = no movement)
Vec	x;
	return inEgg(cell,endN,endC,-len,x);
}

float inEgg ( Vec cell, Seg ends, float len )
{ // just return the distance (-ve len = no movement)
Vec	x;
	return inEgg(cell,ends.A,ends.B,-len,x);
}

int packEgg ( Cell *cell, Cell *child, float strict )
{
Vec	shift;
int	shell,
	level = cell->level,
	type = cell->type,
	sort = cell->sort;
Data	*model = Data::model+cell->model;
float	radius = model->sizes[cell->level]*0.5,
	push   = model->keeps[cell->level]*strict,
	d, dout;
	if (cell->level==Data::depth) return 0;
	if ((int)cell->endN.z == 1234) return 0;	// ends not set yet
	if ((int)cell->endC.z == 1234) return 0;	// ends not set yet
	if (push < -NOISE) shell = 1;			// confine to shell (+/-margin)
		else	   shell = 0;
	d = cell->endN | cell->endC;
	if (d < NOISE) return 0;			// no axis length
	if (cell->ends < 0) { float z = Data::Eratio[sort];
		radius = 0.5*d/z;
	}
	if (shell==0) {
		d = cell->xyz|child->xyz;
		if (d < radius*margin) return 0;	// inside the inscribed sphere
	}
	dout = inEgg(child->xyz,cell->endN,cell->endC,radius*2.0,shift);
	if (shell) {
		if (dout>-0.1 && dout<0.1) return 0;	// close enough to surface
		if (dout>0.0) push = -push;		// switch to pulling in
	} else {
		if (dout < 0.0) return 0;		// inside
	}
	if (dout > 0.0) push *= sqrt(dout);		// go easy on far-away children
	shift *= push;					// push back along resultant
	if (shift.sqr()<NOISE) return 0;
	child->move(shift);
	return 1;
}

/*
sortTube(Cells *cell) 
{
int	i, n = cell->kids;
int	run=0, in=0;
Vec	x,y,z,ave;
	vinit(&ave);
	for (i=0; i<n-1; i++) { Cells *b1,*b2, *p1,*p2;
		b1 = cell->child[i];
		b2 = cell->child[i+1];
		if (b1->link[0]==0) continue;
		if (b2->link[0]==0) continue;
		p1 = b1->link[0];
		p2 = b2->link[0];
		if (b1->id == b2->id-1 && p1->id == p2->id+1) {
			//Pi(b1->id) Pi(p1->id) Pi(b2->id) Pi(p2->id) NL
			vave(b1->xyz,p1->xyz,&x);
			vave(b2->xyz,p2->xyz,&y);
			vsub(x,y,&z);
			vsum(z,&ave);
			run++; in++;
		}
		if (b1->id == b2->id+1 && p1->id == p2->id-1) {
			//Pi(b1->id) Pi(p1->id) Pi(b2->id) Pi(p2->id) NL
			vave(b1->xyz,p1->xyz,&x);
			vave(b2->xyz,p2->xyz,&y);
			vsub(x,y,&z);
			vsum(z,&ave);
			run--; in++;
		}
	}
	vdiv(&ave,(float)in); vnorm(&ave);
	// recheck all against consensus direction (ave)
	for (i=0; i<n-1; i++) { Cells *b1,*b2, *p1,*p2; float d;
		b1 = cell->child[i];
		b2 = cell->child[i+1];
		if (b1->link[0]==0) continue;
		if (b2->link[0]==0) continue;
		p1 = b1->link[0];
		p2 = b2->link[0];
		if (b1->id == b2->id-1 && p1->id == p2->id+1) {
			vave(b1->xyz,p1->xyz,&x);
			vave(b2->xyz,p2->xyz,&y);
			vsub(x,y,&z); vnorm(&z);
			d = vdot(z,ave);
			if (d < -0.2) { // wrong way so swap
				vcopy(b1->xyz,&x); vcopy(b2->xyz,&(b1->xyz)); vcopy(x,&(b2->xyz));
				vcopy(p1->xyz,&x); vcopy(p2->xyz,&(p1->xyz)); vcopy(x,&(p2->xyz));
			}
		}
		if (b1->id == b2->id+1 && p1->id == p2->id-1) {
			vave(b1->xyz,p1->xyz,&x);
			vave(b2->xyz,p2->xyz,&y);
			vsub(x,y,&z); vnorm(&z);
			d = vdot(z,ave);
			if (d < -0.2) { // wrong way so swap
				vcopy(b1->xyz,&x); vcopy(b2->xyz,&(b1->xyz)); vcopy(x,&(b2->xyz));
				vcopy(p1->xyz,&x); vcopy(p2->xyz,&(p1->xyz)); vcopy(x,&(p2->xyz));
			}
		}
	}
}

groupCell (Cells *cell, int level)
{
Vec	centre, shift;
int	i,j,k, m,n, id, type, sort;
float	far = -cell->far/(float)data[3];
float	strict, bias = (float)(cell->level-1);
	if ((int)cell->endN.z==1234 || (int)cell->endC.z==1234) return;
	if (cell->empty) return;
	n = cell->kids;
	if (n==0) return;
	id = cell->id;
	type = cell->type;
	sort = cell->sort;
	if (type<0) type = -type;
	bias = 1.0 - exp(-bias*bias*0.1); // small for top levels to increase selection
	if (level > 0) shifter(cell);
	// readjust centres (except for world) going down
	//
	for (i=0; i<n; i++) {
		groupCell(cell->child[i],level+1);
	}
	//
	// gather lost chlidren back to family coming up
	if (data[0] <= 0) return;
	if (cell->solid > 0) bias = 2.0; else bias = 0.5; // 4x more like to check cthru cell
	if (drand48()*bias > exp(far)) return; // distant cells (small exp(far)) checked less
	// set the model parameters
	model = cell->model;
	if (level==0) model = cell->child[0]->model; // world model type is ambiguous
	m = model*M*N+N;
	moltype = data[m];
	align = data+m+N*1; class = data+m+N*2;
        sizes = data+m+N*3; bumps = data+m+N*4; links = data+m+N*5; chain = data+m+N*6;
	kicks = data+m+N*7; keeps = data+m+N*8; repel = data+m+N*9; bonds = data+m+N*10;
	for (i=0; i<N; i++) keep[i] = 0.01*(float)keeps[i];
	if (rna && type==2 && sort==1) sortTube(cell);
	for (i=0; i<n; i++)
	{ Cells *child = cell->child[i];
	  int	lost = 0; float balance;
	  	if (data[0] < 0) break;  // don't push before setup is done
		if (chain[level] > 0 && chain[level+1] > 0) {   // in a chain of chains (don't push ends in)
			if (child == cell->starts || child == cell->finish) continue;
		}
		strict = 1.0;
		//if (child->type==2 && child->sort==0) continue; // don't push loops
		switch (type) {
			case 0 :	// dummy sphere
				lost = packBall(cell,child,1.0);
			break;
			case 1 :	// simple sphere
				lost = packBall(cell,child,strict);
			break;
			case 2 :	// SSE (tube)
				if (sort && cell->endN.z+cell->endC.z < 1000.0) {
					vave(cell->endN,cell->endC, &(cell->xyz));
				}
//				{ Vec	shift, mid; // bring endpoints in line
//				  float len = 0.1*(float)sizes[level];
//					if (sort==1) len *= 3.0/2.0;
//					vave(cell->endN,cell->endC, &mid);
//					vsub(cell->xyz,mid, &shift);
//					vsum(shift, &(cell->endN));
//					vsum(shift, &(cell->endC));
//					separate(&(cell->endN),&(cell->endC),len,0.1);
//				}
				if (moltype==0) { // protein
					if (sort==0) strict = 0.01; // loose loop
					if (sort==1) strict = 1.00; // tight for alpha
					if (sort==2) strict = 0.10; // looser for beta (to allow bend)
				}
				if (moltype==1) { // nucleic
					strict = 1.0;
				}
				lost = packTube(cell,child,strict);
				if (sort==0) lost = 0; // allow lost children in loops
			break;
			case 3 :	// Ellipsoid
				if (cell->sort < E/2-1 || cell->sort > E/2+1) {
					lost = packEgg(cell,child,1.0);
				} else {	// effectively a sphere
					lost = packBall(cell,child,1.0);
				}
			break;
		}
	}
}
*/
