#define TEST(obj)    if (obj == NULL) { printf("malloc fail for " #obj "\n"); exit(1); }

#include "geom.hpp"
#include "cell.hpp"
#include "data.hpp"
#include <string.h>

//#define MAXIN 100000 //(needed for 3chy/unit3+.run)
#define  MAXIN 20000

void pass1set ( Cell*, int, int, FILE*, int );
void pass2set ( Cell*, int );
void pass3set ( Cell*, int );

Cell    **atom1cell, **atom2cell, **secs2cell;
int      *atom2atom,  *atom2resn,  *resn2atom;

Vec	cent;
char	line[222];
int	secs, atoms, teratoms, endatoms, allatoms;

void models ( FILE *run )
{
int	i,j,k,n, wait;
Cell	*world = new Cell; TEST(world);	// set-up the top cell structure
	secs = teratoms = endatoms = allatoms = 0;
/*
	uid2cell  = (Cells**)malloc(sizeof(Cells*)*MAXIN*8); TEST(uid2cell) // global back-ref.lookup (total)
	atom2cell = (Cells**)alloca(sizeof(Cells*)*MAXIN*6);	// lookup array for all atoms (allatoms)
	secs2cell = (Cells**)alloca(sizeof(Cells*)*MAXIN*1);	// WAS 10000);	// lookup for all SSEs
	atom1cell = (Cells**)alloca(sizeof(Cells*)*MAXIN*1);	// WAS 1000);	// lookup for domain (teratoms)
	atom2resn = (int*)alloca(sizeof(int)*MAXIN*6);	// lookup from global atom id to pdb resid (allatoms)
	atom2atom = (int*)alloca(sizeof(int)*MAXIN*4);	// lookup from local to global atom id (endatoms)
	resn2atom = (int*)alloca(sizeof(int)*MAXIN*1);	// lookup from global pdb resid to atom id
*/
	pass1set(world,0,-1,run,0);
	cent /= (float)allatoms;
	Pi(allatoms) Pv(cent) NL
	pass2set(world,0);
	pass3set(world,0);
//	if (nrelinks) rebond();
	// read cross-link data
//	readlinks(line,pdb); // check to end of file for LINKS
	// set-up the scene data
	world->solid = -1;
	world->live = 999;
/*
	for (i=0; i<SPARE; i++) { // contents allocated in copyCell()
		n = world->kids + i;
		world->child[n]->kids = -1; // flag that no childern exist
		world->child[n]->nbonds = -1; // or bonds to them
		world->child[n]->nlinks = -1; // or links to them
	}
	world->link = uid2cell; // world is never linked so using link to pass uid2cell address
*/
	Cell::world = world;
	printf("Data structures all setup\n");
}

void pass1set ( Cell *cell, int level, int id, FILE *file, int subfile )
// first pass allocates memory, sets values and gets atoms
{
FILE	*pdb;
float	x,y,z, s;
int	i,j,k,m,n, in, sort = 0;
int	nlinks, nbonds, newfile = 0;
Vec	heli, axis, grow, tran, xrot,yrot,zrot;
float	helix, theta, spinx, spiny, spinz;
char	nextfile[111];
int	total, depth = Data::depth;
	helix = theta = spinx = spiny = spinz = 0.0;
	axis.zero(); grow.zero(); tran.zero();
	xrot.zero(); yrot.zero(); zrot.zero();
	cell->atom = cell->ends = 0;
	in = read_line(file,line);
Ps(line) NL
	if (level) cell->xyz = cell->parent->xyz;
	while (line[0]=='#') {	// comments
		Ps(line) NL
		in = read_line(file,line);
	}
	if (line[0]=='M') { // use new MODEL type
		NL Ps(line) NL
		sscanf(line+6,"%d", &model);
		printf("\nSwitching to model %d\n", model);
		moltype = Data::model[m].moltype; // mol.type: 0=PROT, 1=RNA, 2=CHEM
		cell->model = m;
		if (moltype==0) printf("Molecule type = PROTein\n");
		if (moltype==1) printf("Molecule type = RNA\n");
		if (moltype==2) printf("Molecule type = CHEMical\n");
		in = read_line(file,line);
	}
	s = Data::scalein;
	chain = Data::model[model].chain;
	links = Data::model[model].chain;
	nbonds = abs(chain[level]);
	nlinks = links[level];
Pi(nbonds) Pi(nlinks) NL
	if (line[0]=='S' && line[1]=='C') {	// SCALE <x> <y> <z> (applied after all offspring read in)
		Ps(line) NL
		sscanf(line+6,"%f %f %f", &x, &y, &z);
		grow.x = x; grow.y = y; grow.z = z;
		in = read_line(file,line);
	}
	if (line[0]=='H' && line[1]=='E') {	// HELIX <x> <y> <z> <theta> (move x,y, spin helix on Z, move z)
		Ps(line) NL
		sscanf(line+6,"%f %f %f %f", &x, &y, &z, &helix);
		heli.x += x; heli.y += y; heli.z += z;
		helix *= PI/180.0;
		in = read_line(file,line);
	}
	if (line[0]=='T' && line[1]=='W') {	// TWIST <x> <y> <z> <theta> (applied after all offspring read in)
		Ps(line) NL
		axis = cell->xyz;
		sscanf(line+6,"%f %f %f %f", &x, &y, &z, &theta);
		axis.x += x; axis.y += y; axis.z += z;
		theta *= PI/180.0;
		in = read_line(file,line);
	}
	if (line[0]=='S' && line[4]=='X') {	// SPINX <theta> (applied after all offspring read in)
		Ps(line) NL
		sscanf(line+6,"%f", &spinx);
		spinx *= PI/180.0;
		in = read_line(file,line);
	}
	if (line[0]=='S' && line[4]=='Y') {	// SPINY <theta> (applied after all offspring read in)
		Ps(line) NL
		sscanf(line+6,"%f", &spiny);
		spiny *= PI/180.0;
		in = read_line(file,line);
	}
	if (line[0]=='S' && line[4]=='Z') {	// SPINZ <theta> (applied after all offspring read in)
		Ps(line) NL
		sscanf(line+6,"%f", &spinz);
		spinz *= PI/180.0;
		in = read_line(file,line);
	}
	if (line[0]=='S' && line[1]=='H') {	// SHIFT <x> <y> <z> (must be before TWIST in input)
		Ps(line) NL
		sscanf(line+6,"%f %f %f", &x, &y, &z);
		cell->xyz.x += x*s; cell->xyz.y += y*s; cell->xyz.z += z*s;
		in = read_line(file,line);
	}
	if (line[0]=='T' && line[1]=='R') {	// TRANS <x> <y> <z> (applied after offspring read in)
		Ps(line) NL
		sscanf(line+6,"%f %f %f", &x, &y, &z);
		tran.x = x*s; tran.y = y*s; tran.z = z*s;
		in = read_line(file,line);
	}
	if (line[0] == 'I') {	// INPUT[0] <filename> (NB just one space) ("0" resets local atom count)
		Ps(line) NL
		if (line[5]=='0') {
			endatoms = 0;
			Pt(Local ATOM counter set to zero) NL
			strcpy(nextfile,line+7);
		} else {
			strcpy(nextfile,line+6);
		}
		printf("Reading file %s\n", nextfile);
		pdb = fopen(nextfile,"r");
		in = read_line(pdb,line);
		newfile = 1;
		subfile++;
	} else {
		pdb = file;
	}
	cell->endN.zero(); cell->endC.zero();
	cell->endN.z = cell->endC.z = 1234.5;	// marker for unset ends
	cell->cent.x = cell->cent.y = cell->cent.z = -1.0; // -ve dist = unset
	if (line[0] == 'G') {	// GROUP <sort> <members> [<x> <y> <z>  <x> <y> <z> ]
		sscanf(line+6,"%d %d", &sort, &n);
Pi(sort) Pi(n) NL
		if (in > 20) { float x1,y1,z1, x2,y2,z2;
			sscanf(line+11,"%f %f %f  %f %f %f", &x1,&y1,&z1, &x2,&y2,&z2);
			cell->endN.x = x1*s; cell->endN.y = y1*s; cell->endN.z = z1*s;
			cell->endN += cell->xyz;
			cell->endC.x = x2*s; cell->endC.y = y2*s; cell->endC.z = z2*s;
			cell->endC += cell->xyz;
			cell->ends = 1;
		}
		if (level == depth-1) { // SSE level
//			secs2cell[secs] = cell;
			secs++;
		}
	}
	if (line[0] == 'A') // ATOM (PDB format)
	{ int	resn; float a, b;
		sscanf(line+30,"%f %f %f %f %f", &x, &y, &z, &a, &b);
		if (a < NOISE) { // assume unset
			x = drand48()-0.5; y = drand48()-0.5; z = drand48()-0.5;
			if (b > NOISE) { x*=b; y*=b; z*=b; }
		}
		cell->xyz.x  = x*s; cell->xyz.y  = y*s; cell->xyz.z  = z*s;
		cent += cell->xyz;
		sort = n = 0;
		teratoms++; // rezeroed by EoF
		endatoms++; // rezeroed by RELINK/REBOND
		allatoms++; // never rezeroed
		if (allatoms > MAXIN*6) { Pt(Increase MAXIN) Pi(allatoms) NL exit(1); }
		if (endatoms > MAXIN*4) { Pt(Increase MAXIN) Pi(endatoms) Pt(or use INPUT0) NL exit(1); }
		if (teratoms > MAXIN*1) { Pt(Increase MAXIN) Pi(teratoms) NL exit(1); }
		cell->atom = teratoms;		// local number for beta
		cell->btom = endatoms;		// group number for relinks()
		cell->ctom = allatoms;		// number for all
//		atom1cell[teratoms] = cell;	// local to current file
//		atom2cell[allatoms] = cell;	// global over all atoms
//		atom2atom[endatoms] = allatoms; // relink scope to global
		sscanf(line+22,"%d", &resn);
		cell->resn = resn;
		if (resn<0) { Pt(*NB* negative residue numbers are not good: default to 0) NL resn = 0; }
//		resn2atom[resn] = allatoms;
//		atom2resn[allatoms] = resn;
	}
	total = cell->uid; // uid set in spawn()
Pi(total) NL
//	uid2cell[total] = cell;
	if (total > MAXIN*8) { Pt(Increase MAXIN) Pi(total) NL exit(1); }
	if (n==1 && chain[level+1]) { // an only child is not a chain (but might be later)
		printf("*NB* chain of one at level %d\n", level+1);
		// chain[level+1] = 0; 
	}
	cell->kids = cell->cots = n;
	cell->id = id;
	cell->solid = 0;
	cell->far = 9.9;
	for (i=0; i<NDATA; i++) cell->data[i] = 0;
	cell->live = 9; cell->done = cell->lost = cell->bump = cell->busy = cell->empty = 0;
	cell->level = level;
//	cell->type = class[level];	// pick-up the type assigned to the current level
	//if (cell->type==3) sort = E/2;  // default = 1:1 (sphere)
	cell->sort = sort;
	cell->nbonds = nbonds;
/*
	if (nbonds > 0) {
		if (nbonds > 9) { Pt(*NB*) Pi(nbonds) NL exit(1); }
		cell->bond = (Bonds*)malloc(sizeof(Bonds)*nbonds); TEST(cell->bond)
		for (i=0; i<nbonds; i++) {
			cell->bond[i].to = 0;	// bro=0, sis=1
			cell->bond[i].type = 0; // default is unset (0)
			cell->bond[i].next = 0; // default is follow brother
		}
	} else { cell->bond = 0; }
*/
	cell->nlinks = nlinks;
	if (nlinks > 0) {
		cell->link =  new Cell*[nlinks]; TEST(cell->link)
		for (i=0; i<nlinks; i++) cell->link[i] = 0;
	} else { cell->link = 0; }
	for (i=0; i<3; i++) cell->ranks[i] = id;
	cell->junior = cell->senior = 0;
	cell->starts = cell->finish = 0;
	cell->nsis = cell->nbro = 0;
	cell->sis = cell->bro = 0;
	if (n==0) return;
	// make the children
	cell->spawn(n);
	// fill the children
	DO(i,n) pass1set(cell->child[i],level+1,i,pdb,subfile);
	// twists and turns implenented on return path (NB correct input and execution order is critical)
	// HELIX on Z axis: shift by x,y, spin theta on Z then shift z along Z
	if (helix*helix>NOISE && heli.sqr()>NOISE) { float z = heli.z;
		printf("Cell %d at level %d helical rotation by %f rad.\n", cell->id, cell->level, helix);
		moveCell(cell,heli,1);
		heli = cell->xyz; heli.z += 1.0;
		spinCell(cell,cell->xyz,heli,helix);
		heli.x = heli.y = 0.0; heli.z = z;
		moveCell(cell,heli,1);
	}
	// TWIST on general axis (NB axis overwritten in SPIN)
	if (theta*theta > NOISE) {
		printf("Cell %d at level %d rotated by %f rad.\n", cell->id, cell->level, theta);
		spinCell(cell,cell->xyz,axis,theta);
	}
	// SPIN on XYZ axes
	if (spinx*spinx > NOISE) {
		printf("Cell %d at level %d rotated around X by %f rad.\n", cell->id, cell->level, spinx);
		axis = cell->xyz; axis.x += 1.0;
		spinCell(cell,cell->xyz,axis,spinx);
	}
	if (spiny*spiny > NOISE) {
		printf("Cell %d at level %d rotated around Y by %f rad.\n", cell->id, cell->level, spiny);
		axis = cell->xyz; axis.y += 1.0;
		spinCell(cell,cell->xyz,axis,spiny);
	}
	if (spinz*spinz > NOISE) {
		printf("Cell %d at level %d rotated around Z by %f rad.\n", cell->id, cell->level, spinz);
		axis = cell->xyz; axis.z += 1.0;
		spinCell(cell,cell->xyz,axis,spinz);
	}
	// TRANSlate (after spins)
	if (tran.sqr() > NOISE) { float d = tran.len();
		printf("Cell %d at level %d moved by %f\n", cell->id, cell->level, d);
		moveCell(cell,tran,1);
	}
}

void pass2set ( Cell *cell, int level )
// second pass sets links and relationships
{
int	i,j,k, m,n, id, nlinks;
	id = cell->id;
	n = cell->kids;
	m = cell->model;
	moltype = Data::model[m].moltype;
	chain = Data::model[model].chain;
	links = Data::model[model].links;
	nlinks = links[level+1];          // for lower level
	if (moltype==0 && nlinks && cell->type==2) { // SSE
	      	if (cell->sort==1) { // alpha helix
			for (i=0; i<n; i++) { // assume links also allocated at child level
				if (i+3<n)  cell->child[i]->link[0] = cell->child[i+3];
				if (i-4>=0) cell->child[i]->link[1] = cell->child[i-4];
			}
		}
	      	if (cell->sort==2) { // beta sheet
			for (i=0; i<n; i++) { // assume links also allocated at child level
				if (i+2<n)  cell->child[i]->link[0] = cell->child[i+2];
				if (i-2>=0) cell->child[i]->link[1] = cell->child[i-2];
			}
		}
	}
	if (n==0) return;
	// set the last and next children for each child (may be changed by RELINK later)
	m = n-1;
	if (chain[level+1] && n > 1) { // in a chain (NB assumes no mixed models)
		for (i=1; i<m; i++) {
			cell->child[i]->bro = cell->child[i+1];
			cell->child[i]->sis = cell->child[i-1];
		}
		cell->child[0]->bro = cell->child[1];
		cell->child[m]->sis = cell->child[m-1];
		if (chain[level+1] < 0) { // chain is circular so link ends
			cell->child[m]->bro = cell->child[0];
			cell->child[0]->sis = cell->child[m];
		}
		if (chain[level+1] > 0) { // chain is linear so use parent as dummy siblings
			cell->child[m]->bro = cell;
			cell->child[0]->sis = cell;
		}
	}
	if (chain[level+1] && n==1) { // in a 'chain' of one (just dummy siblings)
		cell->child[0]->sis = cell;
		cell->child[m]->bro = cell;
	}
	// set the terminal children (for current cell)
	cell->starts = cell->junior = cell->child[0];
	cell->finish = cell->senior = cell->child[n-1];
	//
	for (i=0; i<n; i++) pass2set(cell->child[i], level+1);
}

void pass3set ( Cell *cell, int level )
// third pass adds links across chain of chains (can only be done after pass2set() )
// sums atoms into higher level centres and sets ends for ellipsoids and tubes
{
int	i,j,k, m,n, id, type, sort;
	id = cell->id;
	n = cell->kids;
	type = abs(cell->type);
	sort = cell->sort;
	m = cell->model;
	moltype = Data::model[m].moltype;
	chain = Data::model[model].chain;
	sizes = Data::model[model].sizes;
	split = Data::model[model].split;
	if (n==0) {
		cell->xyz -= cent;
		return;
	}
	if (n==1) { // one child cannot form a chain (yet)
		if (chain[level+1] != 0) {
			printf("*NB* chain of one at level %d\n", level+1);
			// chain[level+1] = 0;
		}
	}
	if (chain[level] && chain[level+1]>0 && split[level]==0) {
		// cell can be linear or circular but child must be linear and not split (ie forced break)
		if (cell->level == cell->sis->level) {	// set sister of first child as last child of sister
//			cell->starts->bond[1].to = cell->starts->sis = cell->sis->finish;
		}
		if (cell->level == cell->bro->level) {	// set brother of last child as first child of brother
//			cell->finish->bond[0].to = cell->finish->bro = cell->bro->starts;
		}
	}
	for (i=0; i<n; i++) pass3set(cell->child[i], level+1);
	if (level < Data::depth) { // sum children into upper levels
		cell->setWcent();
		if (cell->ends) // ends have been read in (just to set the axis not its ends)
		{ float	size = sizes[level], // for tube and ellipsoids, size = section diameter
		  	axis = cell->endC|cell->endN,
			fn = (float)cell->kids;
		  Vec	x = (cell->endC-cell->endN).norm(); // unit axis in x
/*
			if (type==2) {		// scale SSE tube lengths
				if (moltype==0) { // protein
					if (sort==0) axis = size*sqrt(fn)*loopAXIS;		// loops
					if (sort==1) axis = size*fn*bondCA*alphAXIS/3.8;	// alpha
					if (sort==2) axis = size*fn*bondCA*betaAXIS/3.8;	// beta
				}
				if (moltype==1) { // nucleic
					if (sort==0) axis = 0.2*sqrt(fn);
					// RNA = 2.3 rise/bp --> 0.115 (1/2 for duplex, 1/10 for scale)
					if (sort==1) axis = 0.12*fn;
				}
			}
			if (type==3)	// assign ellipsoid sort
			{ float fn = E/2, nth = 1.0/fn, a =log2(fn+1),
			  	f = (fn-1.0)/(pow(fn,a)*(pow(2.0,a)-1.0)),
				last = 0.0, ratio = axis/size;
			  int	best = -1;
				if (cell->sort==0) {// find preset ellipsoid with best fit to axis/size ratio
					for (i=0; i<E; i++) { float y, fi = (float)i;
						y = nth+f*pow(fi,a);
						if (ratio>last && ratio<y) {
							Pi(i) Pr(last) Pr(y) Pr(ratio) NL
							if (ratio-last < y-ratio) best = i-1; else best = i;
						}
						last = y;
					}
					if (best < 0) { // no fit found
						if (size > axis) {
							best = 0;
							printf("Most oblate ellipsoid used:");
							axis = size/fn;
						} else {
							best = E-1; // for E sorts
							printf("Most prolate ellipsoid used:");
							axis = size*fn;
						}
					} else {	// set axis length to selected ellipsoid
						axis = size*(nth+f*pow((float)best,a));
					}
					Pr(axis) Pr(size) NL
					cell->sort = best;
				} else {
					axis = size*(nth+f*pow((float)cell->sort,a));
				}
			}
*/
			// make symmetric about centre
			x *= axis*0.5;
			cell->endN = cell->xyz - x;
			cell->endC = cell->xyz + x;
			if (moltype==1 && type==2 && sort==1) // for RNA SSE put endC near termini
			{ float dnn, dcn, dnc, dcc;
				dnn = cell->endN|cell->junior->xyz;
				dnc = cell->endC|cell->junior->xyz;
				dcn = cell->endN|cell->senior->xyz;
				dcc = cell->endC|cell->senior->xyz;
				if (dnn+dcn < dnc+dcc) { // swap ends
					cell->endN = cell->xyz + x;
					cell->endC = cell->xyz - x;
				}
			}
			DO(i,n) { Cell *kidi = cell->child[i]; // set internal reference distances
				kidi->cent.x = kidi->xyz|cell->endN;
				kidi->cent.y = kidi->xyz|cell->xyz;
				kidi->cent.z = kidi->xyz|cell->endC;
			}	
		}
	}
}
