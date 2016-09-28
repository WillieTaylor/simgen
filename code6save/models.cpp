#include "util.hpp"
#include "geom.hpp"
#include "cell.hpp"
#include "data.hpp"

//#define MAXIN 100000 //(needed for 3chy/unit3+.run)
#define MAXIN 100000
#define DEEP  10

typedef struct {
	int	grows, helix, twist, shift, spins, trans;
	Vec	grow,  move,  spin,  tran, heli;
	float	turns, theta;
	Seg	axis;
	int	set;
} Moves;

void setCell ( Cell*, int, int, int, int, int );
void setAtom ( Cell*, char* );
void getEnds ( Cell*, char* );
void endCell ( Cell*, Moves* );
int  paired ( Cell*, int, int );
int  moveSet ( char*, Moves* );
void newPath ( char*, char );
int  setModel ( char* );
FILE* getInput ( char*, Moves* );
void pass1set ( FILE* );
void pass2set ( Cell*, int );
void pass3set ( Cell*, int );
void zipDNA ( Cell* );
void setSheet ( FILE*, char );
void setAngle ( Cell*, int );
void readlinks( char*, int );
int  makeLink ( Cell*, Cell*, int, int, int );
int  makeBond ( Cell*, Cell*, int, int, int );
void rebond ();
void rewire ( Cell* );
void setHinge ( char* );

Cell    **atom1cell, **atom2cell, **secs2cell;
int      *atom2atom,  *atom2resn,  *resn2atom;

Vec	cent;
int	secs, atoms, teratoms, endatoms, allatoms;

typedef struct {
        Cell   *at, *to, *scope; // cell source, target, scope of relink
        int     ends;     //  0 = relink, 1 = new start, 2 = new finish
} Relinks;
Relinks *relink;
int     nrelinks;

void models ( FILE *run )
{
int	level;
float	range, size;
Cell	*world = Cell::world;
	world->model = world->level = world->uid = world->id = 0;
	secs = teratoms = endatoms = allatoms = 0;
	model = nrelinks = 0;
        // allocate pairs of cell to collect chain relinks in pass1set
        relink = new Relinks[MAXIN];
	Cell::uid2cell = new Cell*[MAXIN*8]; TEST(Cell::uid2cell) // global back-ref.lookup (total)
	atom2cell = new Cell*[MAXIN*6];	// lookup array for all atoms (allatoms)
	secs2cell = new Cell*[MAXIN*1];	// WAS 10000);	// lookup for all SSEs
	atom1cell = new Cell*[MAXIN*1];	// WAS 1000);	// lookup for domain (teratoms)
	atom2resn = new int[MAXIN*6];	// lookup from global atom id to pdb resid (allatoms)
	atom2atom = new int[MAXIN*4];	// lookup from local to global atom id (endatoms)
	resn2atom = new int[MAXIN*1];	// lookup from global pdb resid to atom id
	Cell::uid2cell[0] = world;
	Cell::total = 0;
	pass1set(run);
	pass2set(world,0);
	pass3set(world,0);
	FOR(i,Cell::total) Cell::uid2cell[i]->done = 0;
	if (nrelinks) rebond(); // remake connections at the atomic level and propagate upwards
	FOR(i,Cell::total) Cell::uid2cell[i]->done = 0;
	level = 0;
	range = 0.1; // must match in main()
	if (Data::frozen) { // a focus has been set
		int uid; Cell *focus;
		if (Data::frozen < 0) { uid = -Data::frozen; } else { uid =  Data::frozen; }
		focus = Cell::uid2cell[uid];
		Data::focus = focus->xyz;
		level = focus->level; // level of focus
		size = 0.5*focus->len;
		printf("Focus = cell:uid %d  level = %d size = %f\n", uid, level, size);
		size *= size;
	} else { // for world (but nothing goes beyond the world anyway)
		size = world->len; size *= size;
	}
	FOR(i,Cell::total+1) // set heat from the squared distance from the focus
	{ Cell	*c = Cell::uid2cell[i];
	  float dd  = c->xyz || Data::focus;
		if (c->level <= level) { // don't hinder levels above the focus
			c->heat = 1.1;	// over 1 = exempt
		} else {
			if (dd < size) { // keep fully active withing focus object
				c->heat = 1.0;
			} else {	// drop-off outside as a Gaussian
				c->heat = exp(-(dd-size)*range/size);
			}
		}
	}
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
	Pv(Data::focus) NL
	Pi(Cell::total) NL
//world->print();
//exit(1);
	FOR(i,world->kids) setAngle(world->child[i],0);
	printf("Data structures all setup\n");
}

void rezero ( Moves *holds ) {
	holds->set = 0;
	holds->grows=holds->helix=holds->twist=holds->shift=holds->spins=holds->trans=0;
	holds->grow.zero(); holds->move.zero();holds->spin.zero();holds->tran.zero();
	holds->turns = holds->theta = 0.0;
}

void pass1set ( FILE *runfile ) {
// read the command <file> to allocate memory, set values and coordinates
/* command syntax
	M = MOVEset = {TRANS, SHIFT, SPIN, HELIX,...} applied to the scope of the GROUP they precede
	G = GROUP   = a group of GROUPs or ATOMs
	A = ATOM    = atomic level coordinate (PDB format)
	T = TER     = end the GROUP (but ends automatically when full (all kids are read)
	B = BONDset = {BOND, LINK, REBOND, REJOIN}  (see below for scope)
  with, lower-case = optional and X* = repeat X, a block has the structure:
	MGATB or mGAtb or mGA*tb or mG*A*t*b or (mG*A*t*b)*
  each time a GROUP is opened, the value of <level> increases and decreases on close
eg: 
-------------------------------------------
level.id (0.0 = world)
1.0	GROUP 0 3
2.0	GROUP 0 6
3.0	ATOM      1  CA  GLY A   1     -23.868 -17.022  -6.339  1.33  1.00
:	:
3.5	ATOM      6  CA  GLY A   6     -17.825 -16.538  -2.703  1.33  1.00
2.1	GROUP 2 5   -14.79 -14.10 -3.15   -4.36 -8.33 -0.83
3.0	ATOM      7  CA  GLY A   7     -14.841 -14.717  -4.363  1.33  2.00
:	:
3.4	ATOM     11  CA  GLY A  11      -4.370  -7.485  -0.937  1.33  2.00
2.2	GROUP 0 6
3.0	ATOM     12  CA  GLY A  12      -1.017  -8.244   0.718  1.33  1.00
:	:
3.5	ATOM     17  CA  GLY A  17      -3.385  -6.721   4.645  1.33  1.00
1.1	GROUP 0 1
2.3	GROUP 0 6
3.0	ATOM      1  CA  GLY A   1     -23.868 -17.022  -6.339  1.33  1.00
:	:
3.5	ATOM      6  CA  GLY A   6     -17.825 -16.538  -2.703  1.33  1.00
-------------------------------------------
FOR(i,10) group[i] = group[0];
Pi(level) Pt(at) FOR(i,4) printf("%3d", at[i]); Pt(uid) FOR(i,4) printf("%3d", group[i]->uid); NL
*/
char	full[222], *line, have, last = 'G'; // the start state
int	level = 0, at[DEEP], setMove, deep, sheet, model = 0;
Cell	*cell, *group[DEEP]; group[0] = Cell::world;
FILE	*files[DEEP], *file;
Moves	moves[DEEP], holds[1];
int	place, baseN, baseC, basepair = 0, print = 0;
	deep = 0;
	files[0] = file = runfile; // current file
	FOR(i,DEEP) { at[i] = moves[i].set = 0; }
	holds->set = 0;
	depth = Data::depth;
	DO // keep reading until all nested GROUPs are complete
	{ int	in = read_line(file,full);
		line = full;
		while (*line == '\t' || *line == ' ') line++;
		have = line[0];
		if (in<0) {	// EoF
			Pt(# EoF) NL
			teratoms = 0;
			if (deep) { // return to upper file stream
				Pt(# Moved back to previous file) NL
				fclose(file);
				file = files[--deep];
				last = 'T';
				continue;
			} else {	// end of runfile (add END)
				strcpy(line,"END"); in = 3;
			}
		}
		if (in<3) {	// junk
			continue;
		}
		if (have=='#') {	// echo comments
			Pt(# Comment) Ps(line+1) NL
			continue;
		}
		Pt(#) Ps(line) NL
		if (strstr(line,"PRINT")) {
			Pt(# Print PDB on exit\n)
			print = 1;
			continue;
		}
		if (strstr(line,"STOP")) {
			printf(" # Exiting\n");
			if (print) sortpdb(1.0);
			exit(1);
		}
		if (line[0]=='R' && line[1]=='E') {	// get a RE[BOND|JOIN|TERM] command
			newPath(line,toupper(line[2]));	// flag =   B    J    T 
			continue;
		}
		if (strstr(line,"TER")) {	// TERminate a group
			teratoms = 0;
			if (last=='T') {
				Pt(# teratoms zeroed by TER command\n)
			}
			continue;
		}
		if (have=='Z') { // zip-up the last (two) DNA chains
			zipDNA(cell->parent);
			continue;
		}
		//
		cell = group[level];	// set current cell for this level
		//
		if (have=='M') {	// switch MODEL type (should occur only between groups)
			model = setModel(line);
			continue;
		}
		if (moveSet(line,holds)) { // get transforms for the next group (store in holds/moves)
			Pt(# Moves) Pi(holds->set) NL
			continue;
		}
		if (have=='I') {     // INPUT[0|*] <filename> [<x> <y> <z>] (NB just one space before file)
			// the sub file must contain a complete group set as EoF = end-of-GROUP. (min = GA)
			files[++deep] = file = getInput(line,holds);
			Pt(#) Pi(deep) NL
			continue;
		}
		if (strstr(line,"SHEET")) { // read (to EOF or END) and set BETA pairs
			Pt(# Beta sheet) NL
			sheet = 1;
			if (in > 6) {
				setSheet(file,toupper(line[6])); // numbering mode specified
			} else {
				setSheet(file,'G'); // default is teratom (group) numbering
			}		
			last = 'S';
			continue;
		}
		if (strstr(line,"BOND") || strstr(line,"LINK")) { // set an atom level bond/link
			readlinks(line,-1);
			continue;
		}
		if (strstr(line,"HINGE")) { //  set a bond between ellipsoid and/or tube ends
			setHinge(line);
			continue;
		}
		if (last=='T') { // close groups up to lowest incomplete level
			DO {
				level--;
				if (level==0) break;
				cell = group[level];
				Pt(# EoG) Pi(level) Pi(cell->uid) Pi(moves[level].set) NL
				endCell(cell,moves+level); // finalises cell settings
				if (strstr(line,"DOUBL") && level==depth-1 ) {
					if (line[6] == '-') break; // NB just one space
				}
				at[level]++;
				if (at[level] < cell->parent->kids) break; // more to fill
			}
			if (level==0) {
				Pt(# all groups complete\n)
			} else { // group[level] = the current cell at each level
				group[level] = cell->parent->child[at[level]];
				cell = group[level];
				Pt(# Next) Pi(cell->uid) Pc(last) Pc(have) NL
			}
		}
		if (strstr(line,"END") || level < 0) {	// END of everything
			Pt(# End\n)
			endCell(group[0],moves);
			Pt(# End of World\n)
			if (print) { // dump the coordinates sorted by GROUP (on internal scale)
				sortpdb(1.0);
				exit(1);
			}
			return;
		}
		if (strstr(line,"DOUBL"))   // switch to nucleic acid base-pairing mode
		{ int	hit, id, n;
			sscanf(line+6,"%d %d", &id, &n);
			hit = paired(cell->parent,id,n);    // scan parent for complementary strand
			if (hit < 0) { // no match (so make new segment)
				basepair = n;
				have = 'G'; 	// follow on as if a new GROUP
			} else { // jump back to the end of the matching complementary strand
				basepair = -n;
				group[level] = cell = cell->parent->child[hit];
				level = cell->level;
				place = at[level];	// remember the current position
				at[level] = hit;	// move to end of matching strand
				last = 'A';
				continue;
			}
		}
		if (strstr(line,"SINGL")) have = 'G';
		if (have=='G')	// GROUP <sort> <members> [<x> <y> <z>  <x> <y> <z> ]
		{ int	n, sort;
			total = Cell::total; // ++ in setCell but Cell::total ++ in spawn()
			sscanf(line+6,"%d %d", &sort, &n);
			if (basepair) // make more slots in parent for basepair SSEs
			{ Cell	*pa = cell->parent; int m; // the parent needs n-1 more slots
				Pt(set BPAIR) Pi(sort) Pi(n) Pi(level) Pi(at[level]) Pi(cell->uid) NL
				m = pa->kids + n-1;	// new family size
				baseN = at[level];	// first strand (current) position in family
				baseC = baseN+n-1;	// final strand position in extended family
				pa->extend(m,sort,0);	// enlarge family and fill from current up
				for (int i=baseN; i<=baseC; i++) { Cell *kidi = pa->child[i];
					setCell(kidi,model,level,sort,i,0);
					kidi->spawn(2);
					FOR(j,2) setCell(kidi->child[j],model,level+1,sort,j,0);
				}
				group[level] = pa->child[baseN];
				pa->child[baseC]->nbro = sort;	// use last SSE to hold id in nbro
				pa->kids = m;		// set new family size
			} else {
				Pt(# set GROUP) Pi(sort) Pi(n) Pi(level) Pi(at[level]) Pi(cell->uid) NL
				setCell(cell,model,level,sort,at[level],n); // fill the current cell
				if (holds->set) moves[level] = *holds; // copy trans on hold to moves
				rezero(holds);
				cell->spawn(n); // make kids and set: id, level, parent and ++Cell::total
				level++;
				if (level==depth) { // set-up cells for all the ATOMs in the GROUP
					FOR(i,n) setCell(cell->child[i],model,level,sort,i,0);
				}
				group[level] = cell->child[0];
				at[level] = 0;
			}
			getEnds(cell,line);	// if more numbers, treat as axis end-points
			sheet = 0;
			last = 'G';
			continue;
		}
		if (have=='A') { int done = 0; // ATOM (PDB format)
	 		Pt(# set ATOM) Pi(level) Pi(at[level]) Pi(cell->uid) NL
			if (basepair) { int s;  // W/C strand = 0/1
				if (basepair > 0) { // on the up strand
					s = 0;
				} else {	  // on the down strand
					s = 1;
				}
				setAtom(cell->child[s],line); // read the atom data
				if (basepair > 0) { // on the up strand
					if (at[level] == baseC) { // end of the strand
						group[depth-1] = cell;
						done = 1;
					} else {
						at[level]++;
					}
				} else {	  // on the down strand
					basepair++;	// increase (-ve) basepair until zero
					if (basepair==0) {	// done (back at first basepair SSE)
						at[depth-1] = place; // reset to current SSE position+1
						group[depth-1] = cell->parent->child[place];
						done = 1;
					} else { 
						at[level]--;
					}
				}
			} else { // not in basepair mode
				if (at[level]==cell->parent->kids-1) done = 1;
				setAtom(cell,line); // read the atom data
				at[level]++;
			}
			if (done) {	// end of family
				level = depth;	// set level to atomic (-1 in basepair mode)
				basepair = 0;
				last = 'T';
			} else {	// move to next
				group[level] = cell->parent->child[at[level]];
				last = 'A';
			}
			continue;
		}
	}
}

int paired ( Cell *cell, int id, int n ) {
// check for an existing complementary strand in the current group
	FOR(i,cell->kids) {
		if (cell->child[i]->nbro == -id) { // found
			Pt(Complementary strand found) NL
			return cell->child[i]->id;
		}
	}
	return -1;
}

int setModel ( char* line ) {
// switch MODEL type (should occur only between groups)
int	model;
	sscanf(line+6,"%d", &model);
	printf("\n # Switching to model %d\n", model);
	moltype = Data::model[model].moltype; // 0=PROT, 1=Nucleic, 2=CHEM, 3=CELL
	subtype = Data::model[model].subtype; // 0=RNA, 1=DNA
	if (moltype==0) printf(" # New molecule type = PROTein\n");
	if (moltype==1) {
		printf(" # New molecule type = Nucleic ");
		if (subtype==1) printf("DNA\n"); else printf("RNA\n");
	}
	if (moltype==2) printf(" # New molecule type = CHEMical\n");
	if (moltype==3) printf(" # New molecule type = CELLs\n");
	return model;
}

int moveSet ( char *line, Moves *moves ) {
Seg	axis;
Vec	grow, move, tran, heli;
float	x,y,z, p,q,r, trans, turns, theta, spinx, spiny, spinz;
float	s = Data::scalein;
	if (line[0]=='S' && line[1]=='C') {	// SCALE <x> <y> <z> (applied after all offspring read in)
		sscanf(line+6,"%f %f %f", &x, &y, &z);
		grow.x = x; grow.y = y; grow.z = z;
		moves->grows = 1; moves->grow = grow;
		moves->set = 1;
		return 1;
	}
	if (line[0]=='H' && line[1]=='E') {	// HELIX <x> <y> <z> <turns> (move x,y, spin about Z, move z)
		sscanf(line+6,"%f %f %f %f", &x, &y, &z, &turns);
		heli.x += x*s; heli.y += y*s; heli.z += z*s;
		turns *= PI/180.0;
		moves->helix = 1; moves->heli = heli; moves->turns = turns;
		moves->set = 2;
		return 1;
	}
	if (line[0]=='T' && line[1]=='W') {	// TWIST <x><y><z> <p><q><r> <theta> (rotate cell around axis)
		sscanf(line+6,"%f%f%f %f%f%f %f", &x,&y,&z, &p,&q,&r, &turns);
		axis.A.x = x*s; axis.A.y = y*s; axis.A.z = z*s;
		axis.B.x = p*s; axis.B.y = q*s; axis.B.z = r*s;
		turns *= PI/180.0;
		moves->twist = 1; moves->axis = axis; moves->turns = turns;
		moves->set = 3;
		return 1;
	}
	if (line[0]=='S' && line[4]=='X') {	// SPINX <theta> (applied after all offspring read in)
		sscanf(line+6,"%f", &spinx);
		spinx *= PI/180.0;
		moves->spins = 1; moves->spin.x = spinx;
		moves->set = 4;
		return 1;
	}
	if (line[0]=='S' && line[4]=='Y') {	// SPINY <theta> (applied after all offspring read in)
		sscanf(line+6,"%f", &spiny);
		spiny *= PI/180.0;
		moves->spins = 2; moves->spin.y = spiny;
		moves->set = 4;
		return 1;
	}
	if (line[0]=='S' && line[4]=='Z') {	// SPINZ <theta> (applied after all offspring read in)
		sscanf(line+6,"%f", &spinz);
		spinz *= PI/180.0;
		moves->spins = 3; moves->spin.z = spinz;
		moves->set = 4;
		return 1;
	}
	if (line[0]=='S' && line[4]=='S') {	// SPINS <theta> (applied after all offspring read in)
		sscanf(line+6,"%f %f %f", &spinx, &spiny, &spinz);
		spinx *= PI/180.0; spiny *= PI/180.0; spinz *= PI/180.0;
		moves->spins = 4; moves->spin = Vec(spinx,spiny,spinz);
		moves->set = 4;
		return 1;
	}
	if (line[0]=='S' && line[2]=='I') {	// SHIFT <x> <y> <z> (must be before TWIST in input)
		sscanf(line+6,"%f %f %f", &x, &y, &z);
		move.x += x*s; move.y += y*s; move.z += z*s;
		moves->shift = 1; moves->move = move;
		moves->set = 5;
		return 1;
	}
	if (line[0]=='T' && line[1]=='R') {	// TRANS <x> <y> <z> (applied after offspring read in)
		sscanf(line+6,"%f %f %f", &x, &y, &z);
		tran.x = x*s; tran.y = y*s; tran.z = z*s;
		moves->trans = 1; moves->tran = tran;
		moves->set = 6;
		return 1;
	}
	return 0; // none found
}

FILE* getInput ( char *line, Moves *moves ) {
// sets up for new call of getGroup() from: INPUT[0|*] <filename> (NB just one space)
char	nextfile[111];
FILE	*file;
char	*more;
float	x, y, z;
int	trans = 0, fat = 6;
	x = y = z = 1.0;
	more = strchr(line+7,' '); // starts at any text after the filename
	if (more) {	// read the position (saves having separate TRANS line)
		sscanf(more+1,"%f %f %f", &x, &y, &z);
		*more = (char)0;	// stop x,y,z being part of filename
		trans = 1;
	}
	if (line[5]=='0') { // "0" resets local atom counts
		teratoms = endatoms = 0;
		Pt(# Local ATOM counter set to zero) NL
	}
	if (line[5]=='*') { // "*" adds a random scatter (default = isotropic)
		x = x*(randf()-0.5);
		y = y*(randf()-0.5);
		z = z*(randf()-0.5);
		trans = 1;
	}
	if (trans) {
		moves->tran.x = x; moves->tran.y = y; moves->tran.z = z;
		moves->set = moves->trans = 1;
	}
	if (line[5] != ' ') fat++;
	strcpy(nextfile,line+fat);
	printf(" # Reading file %s\n", nextfile);
	file = fopen(nextfile,"r");
	return file;
}

void setCell ( Cell *cell, int model, int level, int sort, int id, int n ) {
// fill the <cell> with data (except don't know number of children until GROUP command)
int	nbonds =  Data::model[model].chain[level],
	nlinks =  Data::model[model].links[level],
	type   =  Data::model[model].shape[level];
float	size   =  Data::model[model].sizes[level];
	total++;
	if (total > MAXIN*8) { Pt(Increase MAXIN) Pi(total) NL exit(1); }
	if (n==1 && Data::model[model].chain[level+1]) { // an only child is not a chain
		printf("*NB* chain of one at level %d\n", level+1);
	}
	if (level == depth-1) { // SSE level
		secs2cell[secs] = cell;
		secs++;
	}
	cell->ends = 0;
	cell->endN.zero(); cell->endC.zero();
	cell->endN.z = cell->endC.z = 1234.5;	// marker for unset ends
	cell->cent.x = cell->cent.y = cell->cent.z = -1.0; // -ve dist = unset
	cell->kids = cell->cots = n;
	cell->model = model;
	cell->solid = 0;
	cell->far = 9.9;
	cell->len = size;	// default length is the diameter
	cell->heat = 0.0;
	FOR(i,NDATA) { cell->idata[i] = 0; cell->fdata[i] = 0.0; cell->vdata[i].zero(); }
	cell->live = 9;
	cell->atom = cell->btom = cell->ctom = 0;
	cell->done = cell->lost = cell->bump = cell->busy = cell->empty = 0;
	cell->hit = 0;
	cell->type = type;	// pick-up the type assigned to the current level
	cell->sort = sort;
	if (cell->type==3 && sort==0) sort = E/2;  // default = 1:1 (sphere)
	// set the bonds
	cell->nbonds = nbonds;
	if (nbonds < 0) nbonds = -nbonds;	// -ve flags circular topology
	// sis and bro are shorthand for bond[0].to and bond[1].to (and exist even if no bonds)
	cell->sis = cell->bro = cell->parent;
	if (nbonds > 0) { int i;
		if (nbonds > 9) { Pt(*NB*) Pi(nbonds) NL exit(1); }
		cell->bond = new Bonds[nbonds+1]; TEST(cell->bond) // +1 to cover sis=0 + bro=1
		for (i=0; i<=nbonds; i++) {	// NB: <=
			cell->bond[i].to = 0;	// unset
			cell->bond[i].type = 0; // default is thickness 0 (line) 
			cell->bond[i].link = 0; // default connection is N--C (-1 = don't link)
			cell->bond[i].next = i; // default is go to same
		}
		cell->bond[0].to = cell->bond[1].to = cell->parent; 
	} else { cell->bond = 0; }
	// set the links
	cell->nlinks = nlinks;
	if (nlinks > 0) { int i;
		cell->link =  new Bonds[nlinks]; TEST(cell->link)
		for (i=0; i<nlinks; i++) {
			cell->link[i].to = 0;	// unset
			cell->link[i].type = 0; // default is thickness 0 (line) 
			cell->link[i].next = 5; // default is break when over 5% over extended
		}
	} else { cell->link = 0; }
	FOR(i,3) cell->ranks[i] = id;
	cell->junior = cell->senior = 0;
	cell->starts = cell->finish = 0;
	cell->nsis = cell->nbro = 0;
}

void setAtom ( Cell *cell, char *line ) {
// read the coordinates and set atom counters
int	resn; float x,y,z, a,b, s = Data::scalein;
	sscanf(line+30,"%f %f %f %f %f", &x, &y, &z, &a, &b);
	if (a < NOISE) { // assume unset
		x = drand48()-0.5; y = drand48()-0.5; z = drand48()-0.5;
		if (b > NOISE) { x*=b; y*=b; z*=b; }
	}
	cell->xyz.x = x*s; cell->xyz.y = y*s; cell->xyz.z = z*s;
	teratoms++; // rezeroed by EoF or TER 
	endatoms++; // rezeroed by REJOIN/REBOND
	allatoms++; // never rezeroed
	if (allatoms > MAXIN*6) { Pt(Increase MAXIN) Pi(allatoms) NL exit(1); }
	if (endatoms > MAXIN*4) { Pt(Increase MAXIN) Pi(endatoms) Pt(or use INPUT0) NL exit(1); }
	if (teratoms > MAXIN*1) { Pt(Increase MAXIN) Pi(teratoms) NL exit(1); }
	cell->atom = teratoms;		// local number for beta
	cell->btom = endatoms;		// group number for relinks()
	cell->ctom = allatoms;		// number for all
	atom1cell[teratoms] = cell;	// local to current file
	atom2cell[allatoms] = cell;	// global over all atoms
	atom2atom[endatoms] = allatoms; // relink scope to global
	sscanf(line+22,"%d", &resn);
	cell->resn = resn;
	if (resn<0) { Pt(*NB* negative residue numbers are not good: default to 0) NL resn = 0; }
	resn2atom[resn] = allatoms;
	atom2resn[allatoms] = resn;
	// useful data for driver()
	cell->idata[0] = acid2num(line+17);
	strncpy(cell->cdata,line+17,3); cell->cdata[3] = (char)0;
	if (line[21] != ' ') { // set chain id in cdata[4] (and flag it in idata[4])
		cell->cdata[4] = line[21]; // chain ID
		cell->idata[4] = 1;
	}
	cell->fdata[0] = a;
	cell->fdata[1] = b;
}

void endCell ( Cell *cell, Moves *moves ) {
int	n = cell->kids, sort = cell->sort, m = cell->model;
Data	*param = Data::model+m;
int	moltype = param->moltype,
	subtype = param->subtype;
	Pt(# End of Cell) Pi(cell->level) Pi(cell->id) Pi(cell->uid) Pi(cell->kids) Pi(cell->ends) NL
	// set cell position to average of children
	cent.zero();
	FOR(i,n) cent += cell->child[i]->xyz;
	cell->xyz = cent/(float)n;
	if (cell->ends) { // shift axis poles (vdata set in getEnds())
		cell->endN = cell->xyz - cell->vdata[0];
		cell->endC = cell->xyz + cell->vdata[0];
	} else {
		if (moltype==0 && cell->level == depth-1) { // Protein SSE with no ends
			if (n==2) { // use termini
				cell->endN = cell->child[0]->xyz; 
				cell->endC = cell->child[1]->xyz; 
			}
			if (n>2) { // use average of first/last pair (good for beta and loops)
				cell->endN = cell->child[0]->xyz & cell->child[1]->xyz; 
				cell->endC = cell->child[n-1]->xyz & cell->child[n-2]->xyz; 
			}
			if (n>3 && cell->sort==1) { // alpha
				cell->endN = cell->endN & cell->child[3]->xyz; 
				cell->endC = cell->endC & cell->child[n-3]->xyz; 
			}
			cell->ends = 1;
		}
	}
	// apply transformations to current cell and contents
	if (moves->set) { float x,y,z; Vec vec; Seg seg;
		if (moves->helix) { // HELIX rotate about <heli> and move by <heli>
			moves->helix = 0;
			vec = moves->heli; vec.z = 0;
			cell->move(vec);		// shift off axis
			seg.B.x = seg.B.y = 0; seg.B.z = 1;	// Z axis line
			cell->spin(seg,moves->turns);	// spin about Z
			vec.x = vec.y = 0; vec.z = moves->heli.z;
			cell->move(vec);		// shift along Z
			Pt(# Helical twist along Z) Pv(moves->heli) Pr(moves->turns) NL
		}
		if (moves->twist) { // TWIST rotate about <axis> by <twist>
			moves->twist = 0;
			cell->spin(moves->axis,moves->turns);
			Pt(# Twist about) Pl(moves->axis) Pr(moves->turns) NL
		}
		if (moves->spins) { // SPINS about XYZ axes
			moves->spins = 0;
			x = moves->spin.x;
			if (x*x > NOISE) cell->spin(Vec(1,0,0),x);
			y = moves->spin.y;
			if (y*y > NOISE) cell->spin(Vec(0,1,0),y);
			z = moves->spin.z;
			if (z*z > NOISE) cell->spin(Vec(0,0,1),z);
			Pt(# Rotated about XYZ by) Pv(moves->spin) NL
			
		}
		if (moves->trans) { // TRANSlate
			moves->trans = 0;
			Pt(#) Pv(moves->tran) NL
			cell->move(moves->tran);
			Pt(# Moved to) Pv(cell->xyz) NL
		}
		moves->set = 0;
	}
}

void getEnds ( Cell *cell, char *line ) {
// read axis (if given)
float	x1,y1,z1, x2,y2,z2; Vec v;
	if (strlen(line) < 20) return;
	sscanf(line+12,"%f %f %f  %f %f %f", &x1,&y1,&z1, &x2,&y2,&z2);
	v = Vec(x2-x1, y2-y1, z2-z1);
	v *= 0.5*Data::scalein;
	cell->vdata[0] = v; // temp store then N-C split across centre in endCell()
	cell->ends = 1;
	if (cell->sort < 0) {// -ve = flag to use the given axis length
		cell->sort = -cell->sort; cell->ends = -1;
	}
	Pt(# Axis ends read) Pv(v) Pi(cell->ends) NL
}

void newPath ( char *line, char mode ) {
// get the rewiring moves for a new bond path (mode: B=bond, J=join, T=term)
int	at, to; Cell *cat, *cto;
	sscanf(line+13,"%d %d", &at, &to);
	if (toupper(line[7])=='P') { // use pdbid
		cat = atom2cell[resn2atom[at]];
		cto = atom2cell[resn2atom[to]];
	}
	if (toupper(line[7])=='G') { // use (group) id
		cat = atom1cell[at];
		cto = atom1cell[to];
	}
	if (toupper(line[7])=='L') { // local numbering (seq. from: start, INPUT0 or REBOND)
		cat = atom2cell[atom2atom[at]]; // atom2atom[] gives the global (allatom) id 
		cto = atom2cell[atom2atom[to]];
	}
	if (toupper(line[7])=='A') { // atom level numbering
		cat = atom2cell[at];
		cto = atom2cell[to];
	}
	if (toupper(line[7])=='U') { // unique (uid) numbering
		cat = Cell::uid2cell[at];
		cto = Cell::uid2cell[to];
	}
	// in T mode, cat = new N-terminus, cto = new C-terminus
	relink[nrelinks].at = cat;
	relink[nrelinks].to = cto;
	//Pt(New path) Pc(mode) Pi(cat->uid) Pi(cto->uid) Pi(nrelinks) NL
	if (mode=='T') relink[nrelinks].ends = 0;
	if (mode=='B') relink[nrelinks].ends = 1;
	if (mode=='J') relink[nrelinks].ends = 2;
	nrelinks++;
	teratoms = endatoms = 0;
}

#define SNAP 50 // snap link if stretched beyond 50%

void pass2set ( Cell *cell, int level )
// second pass sets chains, links and relationships
{
int	kidbonds;
int	i,j,k, m,n, id, nlinks;
	id = cell->id;
	n = cell->kids;
	m = cell->model;
	moltype = Data::model[m].moltype;
	subtype = Data::model[m].subtype;
	chain = Data::model[m].chain;
	links = Data::model[m].links;
	nlinks = links[level+1];          // for lower level
	if (nlinks && moltype==0 && cell->type==2) { // SSE (link->next is used as %x10 break point)
		if (nlinks < 2) { Pi(nlinks) Pt(SSE children need at least 2 links\n)  exit(1); }
		for (i=0; i<n; i++) if (cell->child[i]->nlinks < 2) cell->child[i]->nlinks = 2; // may be set
	      	if (cell->sort==1) { // alpha helix
			FOR(i,n) // foreach atom (i) in the helix set i+3, i+4 links
			{ Cell	*kidi = cell->child[i];
				if (nlinks < 3) Pt(NB: links overwritten by alpha\n)
				if (kidi->link->to != 0) { // something already linked (so bump-up existing)
					for (j=nlinks-1; j>1; j--) {
						kidi->link[j] = kidi->link[j-2];
					}
				}
				if (i+3<n)  { Bonds *link = kidi->link+0;
					link->to = cell->child[i+3];
					link->type = 10;
					link->next = SNAP;
				}
				if (i-4>=0) { Bonds *link = kidi->link+1;
					link->to = cell->child[i-4];
					link->type = 11;
					link->next = SNAP;
				}
			}
		}
	      	if (cell->sort==2) { // beta sheet
			FOR(i,n) // foreach atom (i) in the strand set i+2, i-2 links
			{ Cell	*kidi = cell->child[i];
				if (kidi->link->to != 0) { // something already linked (so bump-up)
					if (nlinks < 5) Pt(NB: links overwritten by beta\n)
					for (j=nlinks-1; j>3; j--) {
						kidi->link[j] = kidi->link[j-4];
					}
				}
				if (i+2<n)  { Bonds *link = kidi->link+0;
					link->to = cell->child[i+2];
					link->type = 12;
					link->next = SNAP;
				}
				if (i-2>=0) { Bonds *link = kidi->link+1;
					link->to = cell->child[i-2];
					link->type = 13;
					link->next = SNAP;
				}
			}
		}
	}
	if (n==0) return;
	// set the last and next children for each child (may be changed by REJOIN later)
	if (cell->child[0]->bond==0) kidbonds = 0; else kidbonds = 1;
	m = n-1;		// m = index of last child[0...m]
	if (chain[level+1] && n > 1) { // in a chain (NB assumes no mixed models)
		FOR(i,n)
		{ Cell	*kidi = cell->child[i];
		  int	modi = kidi->model;
			if (i<m) { int modj = cell->child[i+1]->model;
				if (modi==modj) kidi->bro = cell->child[i+1];
			}
			if (i>0) { int modh = cell->child[i-1]->model;
				if (modh==modi) kidi->sis = cell->child[i-1];
			}
			if (!kidbonds) continue;
			kidi->bond[0].to = kidi->sis;
			kidi->bond[1].to = kidi->bro;
			kidi->nbro = 1; // only count brothers (always just one sister)
		}
		if (chain[level+1] < 0) { // chain is circular so link ends
			cell->child[m]->bro = cell->child[0];
			cell->child[0]->sis = cell->child[m];
			if (kidbonds) {
				cell->child[m]->bond[1].to = cell->child[0];
				cell->child[m]->bond[1].next = 1;
				cell->child[0]->bond[0].to = cell->child[m];
				cell->child[0]->bond[0].next = 0;
			}
		}
		if (chain[level+1] > 0) { // chain is linear so use parent as dummy siblings
			cell->child[m]->bro = cell;
			cell->child[0]->sis = cell;
			if (kidbonds) {
				cell->child[m-1]->bond[1].next = 0; // cyclise the bond path (only at
				cell->child[ 1 ]->bond[0].next = 1; // the ends of an unbranched chain)
			}
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
	if (moltype==1 && level==depth-2 && chain[level+1] && chain[level+2]) // in a NA seg.
	{ int	seg[4]; // for termini of double-stranded (ds) segments
		FOR(i,n) // loop over SSE segments to set internal ss and ds bonds (sort: ss=0, ds=1)
		{ Cell	*kidi = cell->child[i]; // kidi = a base-paired rung or loop (SSE level)
			if (kidi->sort==0) { // single-stranded loop segment
				// ??? DO1(j,kidi->kids-1) { // set internal loop bonds
				FOR(j,kidi->kids) { Cell* cj = kidi->child[j]; // set internal loop bonds
					if (j) {
						cj->sis = cj->bond[0].to = kidi->child[j-1];
					}
					if (j<kidi->kids-1) {
						cj->bro = cj->bond[1].to = kidi->child[j+1];
					}
				}
				kidi->junior = kidi->starts = kidi->child[0];
				kidi->senior = kidi->finish = kidi->child[kidi->kids-1];
				kidi->endN = kidi->junior->xyz; kidi->endC = kidi->senior->xyz;
				kidi->ends = 1;
			} else		// double-stranded base-paired segment
			{ Cell	*W = kidi->child[0], *C = kidi->child[1]; // Watson/Crick strands
				if (i && cell->child[i-1]->sort!=0) {
					W->sis = W->bond[0].to = cell->child[i-1]->child[0];
					C->bro = C->bond[1].to = cell->child[i-1]->child[1];
				}
				if (i<n-1 && cell->child[i+1]->sort!=0) {
					C->sis = C->bond[0].to = cell->child[i+1]->child[1];
					W->bro = W->bond[1].to = cell->child[i+1]->child[0];
				}
				kidi->endN = W->xyz; kidi->endC = C->xyz;
				kidi->junior = kidi->starts = W;
				kidi->senior = kidi->finish = C;
				kidi->ends = 1;
			}
		}
		FOR(i,4) seg[i] = -1;
		DO1(i,n-1) // loop over the SSE termini and set bonds between ds and ss segments
		{ Cell *kidi = cell->child[i-1], *kidj = cell->child[i];
		  int	si = abs(kidi->sort), sj = abs(kidj->sort);
			if (si>0 && sj>0) { // pair..pair (in a ds seg)
				if (seg[0]<0) { // ds seg start
					seg[0] = kidi->child[0]->ctom;
					seg[1] = kidi->child[1]->ctom;
				}
				seg[2] = kidj->child[0]->ctom;
				seg[3] = kidj->child[1]->ctom;
				continue;
			}
			if (seg[0]<0) continue; // no ds segs yet
			FOR(k,n) // check all seg ends against current ds ends (0=Wn, 1=Cc, 2=Wc, 3=Cn)
			{ Cell	*kid = cell->child[k],
				*jun = kid->junior, *sen = kid->senior, *sej;
			  int	nk = jun->ctom, ck = sen->ctom;
				if (kid->sort != 0) continue;
				FOR(j,4) if (nk==seg[j]+1) {
					sej = atom2cell[seg[j]];
					Pt(N) Pi(k) Pi(j) Pi(seg[j]) Pi(nk) Pi(jun->resn) Pi(sej->resn) NL
					jun->sis = jun->bond[0].to = sej;
					sej->bro = sej->bond[1].to = jun;
					Pi(jun->sis->resn) Pi(sej->bro->resn) NL
				}
				FOR(j,4) if (nk==seg[j]-1) {
					sej = atom2cell[seg[j]];
					Pt(N) Pi(k) Pi(j) Pi(seg[j]) Pi(nk) Pi(sej->resn) Pi(jun->resn) NL
					sej->sis = sej->bond[0].to = jun;
					jun->bro = jun->bond[1].to = sej;
					Pi(sej->sis->resn) Pi(jun->bro->resn) NL
				}
				FOR(j,4) if (ck==seg[j]+1) {
					sej = atom2cell[seg[j]];
					Pt(C) Pi(k) Pi(j) Pi(seg[j]) Pi(ck) Pi(sen->resn) Pi(sej->resn) NL
					sen->sis = sen->bond[0].to = sej;
					sej->bro = sej->bond[1].to = sen;
					Pi(sen->sis->resn) Pi(sej->bro->resn) NL
				}
				FOR(j,4) if (ck==seg[j]-1) {
					sej = atom2cell[seg[j]];
					Pt(C) Pi(k) Pi(j) Pi(seg[j]) Pi(ck) Pi(sej->resn) Pi(sen->resn) NL
					sej->sis = sej->bond[0].to = sen;
					sen->bro = sen->bond[1].to = sej;
					Pi(sej->sis->resn) Pi(sen->bro->resn) NL
				}
			}
			FOR(j,4) seg[j] = -1;
		}
	} else { // continue down
		FOR(i,n) pass2set(cell->child[i], level+1);
	}
}

void pass3set ( Cell *cell, int level )
// third pass adds links across chain of chains (can only be done after pass2set() )
// sums atoms into higher level centres and sets ends for ellipsoids and tubes
{
float	bondCA = Data::bondCA;
int	i,j,k, dna,
	m = cell->model,
	n = cell->kids,
	id = cell->id,
	type = abs(cell->type),
	sort = cell->sort;
Data	*param = Data::model+m;
	moltype = param->moltype;
	subtype = param->subtype;
	chain = param->chain;
	sizes = param->sizes;
	split = param->split;
	dna = moltype;//*subtype; // 1 = true
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
	if (chain[level] && chain[level+1]>0 && split[level]==0 && dna==0) { // join cousins across parents
		// cell can be linear or circular but child must be linear and not split (ie forced break)
		if (cell->level == cell->sis->level) {	// set sister of first child as last child of sister
			cell->starts->sis = cell->sis->finish;
			if (cell->starts->bond) cell->starts->bond[0].to = cell->sis->finish;
		}
		if (cell->level == cell->bro->level) {	// set brother of last child as first child of brother
			cell->finish->bro = cell->bro->starts;
			if (cell->finish->bond) cell->finish->bond[1].to = cell->bro->starts;
		}
	}
	//
	for (i=0; i<n; i++) pass3set(cell->child[i], level+1);
	//
	if (level < Data::depth)
	{ float	size = sizes[level],	// for tube and ellipsoids, size = section diameter
		fn = (float)cell->kids,
		axis = 1.0;
	  Vec	x;
		cell->setWcent();	// sum (weighted) children into upper levels
		if (cell->ends) {	// ends have been read in
			axis = cell->endC | cell->endN;	// set length for ellipsoid (unless sort>0 or ends<0)
			x = cell->endC - cell->endN;		// set axis direction
		} else {
			if (cell->kids > 1) {	// take the direction between termini
				x = cell->senior->xyz - cell->junior->xyz;
			} else {
				x.set_rand(); 	// pick a random axis direction
			}
		}
		x.setVec();		// set to unit length
		if (type==1) {		// for sphere, lenght = diameter
			axis = cell->len;
		}
		if (type==2)		// scale SSE tube lengths
		{ float	sf = fn-1.0;
			if (cell->ends > 0) { // use ideal axis length (otherwise keep the given length)
				if (moltype==0) { // protein
					if (sort==0) axis = loopAXIS*size*sqrt(sf); // loops
					if (sort==1) axis = 0.7*alphAXIS*sf*bondCA/3.8; // alpha
					if (sort==2) axis = betaAXIS*sf*bondCA/3.8; // beta
				}
				if (moltype==1) {
					if (subtype==0 && level==depth-2) { // RNA tube
						if (sort==0) axis = 0.2*sqrt(fn);
						// RNA = 2.3 rise/bp --> 0.115 (1/2 for duplex, 1/10 for scale)
						if (sort==1) axis = 0.12*fn;
					}
					if (level==depth-1) { // DNA/RNA (rungs of basepaired ladder)
						axis *= 0.7; // keep ends inside Ps
					} 
				}
			}
		}
		if (type==3)	// assign ellipsoid sort
		{ float last = 0.0, ratio = axis/size;
		  int	best = -1;
			if (cell->sort==0) {// find preset ellipsoid with best fit to axis/size ratio
				for (i=1; i<E; i++) { float y, fi = (float)i;
					y = Data::Eratio[i];
					if (ratio>last && ratio<y) {
						Pi(i) Pr(last) Pr(y) Pr(ratio) NL
						if (ratio-last < y-ratio) best = i-1; else best = i;
					}
					last = y;
				}
				if (best < 0) { // no fit found
					if (size > axis) {
						best = 1;
						printf("Most oblate ellipsoid used:");
					} else {
						best = E-1; // for E sorts
						printf("Most prolate ellipsoid used:");
					}
				} 	// set axis length to selected ellipsoid
				axis = size*Data::Eratio[best];
				Pi(best) Pr(axis) Pr(size) NL
				cell->sort = best;
			} else { // sort given. If ends<0 keep the given axis length otherwise use the ideal length
				if (cell->ends >= 0) axis = size*Data::Eratio[cell->sort];
			}
		}
		// make symmetric about centre
		x *= axis*0.5;
		cell->endN = cell->xyz - x;
		cell->endC = cell->xyz + x;
		if (moltype==1 && subtype==0 && type==2 && sort==1) // for RNA SSE put endC near termini
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
		FOR(i,n) { Cell *kidi = cell->child[i]; // set internal reference distances
			kidi->cent.x = kidi->xyz|cell->endN;
			kidi->cent.y = kidi->xyz|cell->xyz;
			kidi->cent.z = kidi->xyz|cell->endC;
		}	
		//cell->len = axis; /////////////////////////////////
	}
}

void zipDNA( Cell *cell ) {
int	n = cell->kids;
	FOR(i,n/2)
	{ Cell	*w = cell->child[i],
		*c = cell->child[n-i-1];
	  float d = w->xyz|c->xyz;
		w->link[0].to = c;
		c->link[0].to = w;
	}
}

void setSheet ( FILE *beta, char mode )
{ // set the beta sheet links (link[].next = c has weight c/10)
char	line[222];
int	nlinks;
	while (1) { int i, io, a,b,c; Cell *from, *to;
		io = read_line(beta,line);
		if (io < 0 ) break;
		if (io < 10) continue;
		if (line[0] == '#') continue;
		// BETA 123 456  5
		sscanf(line+5,"%d %d %d", &a, &b, &c); // c = strength
		if (c>5) { Pt(#) Ps(line) NL } // echo strong ones
		if (mode=='G') {	// set links using group ids
			from = atom1cell[a]; to = atom1cell[b];
		}
		if (mode=='A') {	// set links using atom ids
			from = atom2cell[a]; to = atom2cell[b];
		}
		if (mode=='U') {	// set links using unique ids
			from = Cell::uid2cell[a]; to = Cell::uid2cell[b];
		}
		if (mode=='P') { // use pdbid
			from = atom2cell[resn2atom[a]];
			to = atom2cell[resn2atom[b]];
		}
		if (mode=='L') { // local numbering (seq. from: start, INPUT0 or REBOND)
			from = atom2cell[atom2atom[a]]; // atom2atom[] gives the global (allatom) id 
			to = atom2cell[atom2atom[b]];
		}
		//if (from->parent->sort<2 || to->parent->sort<2) continue; // SHEET atoms not in beta
		nlinks = Data::model[from->model].links[from->level];
		for (i=2; from->link && i<nlinks; i++) { // 0,1 used for +/-2 in strand
			if (from->link[i].to == to) break; // already set
			if (from->link[i].to == 0) {
				Pt(# set BETA) Pi(i) Pi(from->resn) Pi(to->resn) NL
				from->link[i].to = to;
				from->link[i].type = 15;
				from->link[i].next = -c; // -ve = not brittle (used as wt in fixSheet)
				from->nlinks = i+1;
				break;
			}
		}
		from = from->parent; to = to->parent;	// set links at strand (midpoint) level
		nlinks = Data::model[from->model].links[from->level];
		if (from->sort<2 || to->sort<2) continue; // some SHEET2 atoms may be in loops
		for (i=0; from->link && i<nlinks; i++) {
			if (from->link[i].to == to ) break;  // already set
			if (from->link[i].to == 0) {
				from->link[i].to = to;
				from->link[i].type = 9;
				from->link[i].next = 1;
				from->nlinks = i+1;
				break;
			}
		}
	}
}

void setAngle ( Cell *cell, int update )
{
// alpha theta, tau = 90.5, 50.4
// alpha theta, tau =  1.6, -0.9
// betas theta, tau =  2.2,  2.7
//
// alpha d 0-2,3 = 5.4, 5.1, 6.2
// betas d 0-2,3 = 6.7, 9.9,13.5
//
float	theta, tau1, tau2;
float	dis12, dis1, dis2;
int	i, n, id, level = cell->level;
float	f, g, fix = 0.1;
Vec	alphA, alphD, betaA, betaD;
float	s = Data::bondCA/3.8; // scale internal/real
	alphA.x = alphA.z = -0.9; alphA.y =  1.6;
	betaA.x = betaA.z =  2.7; betaA.y =  2.2;
	alphD.x = alphD.z = 5.1*s; alphD.y =  6.2*s;
	betaD.x = betaD.z = 9.9*s; betaD.y = 13.5*s;
	id = cell->id;
	n = cell->kids;
	theta = tau1 = tau2 = 999.9;
	dis12 = dis1 = dis2 = 999.9;
	if (cell->sis && cell->bro && cell->parent->kids>3) { // connected at both ends in a chain
		theta = angle(cell->sis->xyz, cell->xyz, cell->bro->xyz);
		if (theta > NOISE) {
			if (cell->sis->sis) {
				if (angle(cell->sis->sis->xyz, cell->sis->xyz, cell->xyz) > NOISE) {
					tau1 = torsion(cell->sis->sis->xyz,cell->sis->xyz,cell->xyz,cell->bro->xyz);
				}
				dis1 = vdif(cell->sis->sis->xyz,cell->bro->xyz); // i-2..i+1
			}
/*
			if (cell->bro->bro) {
				if (angle(cell->xyz, cell->bro->xyz, cell->bro->bro->xyz) > NOISE) {
					tau2 = torsion(cell->bro->bro->xyz,cell->bro->xyz,cell->xyz,cell->sis->xyz);
				}
				dis2 = vdif(cell->bro->bro->xyz,cell->sis->xyz); // i+2..i-1
			}
*/
		}
		if (dis1<999.0 && dis2<999.0 && cell->parent->kids>4) {
			dis12 = vdif(cell->bro->bro->xyz,cell->sis->sis->xyz); // i+2..i-2
		}
		if (update > 0) {
			f = fix; g = 1.0-f;
			cell->geom.x = g*cell->geom.x + f*tau1;
			cell->geom.y = g*cell->geom.y + f*theta;
			cell->geom.z = g*cell->geom.z + f*tau2;
			//f = fix*fix; g = 1.0-f;
			cell->dist.x = g*cell->dist.x + f*dis1;
			cell->dist.y = g*cell->dist.y + f*dis12;
			cell->dist.z = g*cell->dist.z + f*dis2;
		} else {
			cell->geom.x = tau1; cell->geom.y = theta; cell->geom.z = tau2;
			cell->dist.x = dis1; cell->dist.y = dis12; cell->dist.z = dis2;
		}
		if (cell->sis->type==2 && cell->type==2 && cell->bro->type==2) { // in SSE (not end)
			if (cell->sis->sort==1 && cell->sort==1 && cell->bro->sort==1) { // in alpha
				vave(cell->geom, alphA, &(cell->geom));
				vave(cell->dist, alphD, &(cell->dist));
			}
			if (cell->sis->sort==2 && cell->sort==2 && cell->bro->sort==2) { // in beta
				vave(cell->geom, betaA, &(cell->geom));
				vave(cell->dist, betaD, &(cell->dist));
			}
		}
	}
	if (level && cell->bond && update<1) // local bond lengths not updated
	{ int	m = cell->model;
	  Data	*param = Data::model+m;
	  float bond = param->sizes[level]+param->bonds[level];
	  Cell	*sis = cell->sis, *bro = cell->bro;
	  int	freesis, freebro; // -ve dist = free bond (set by REJOIN in rebond())
	  	freesis = freebro = 0;
	  	if (cell->prox.x < -NOISE) freesis = 1;
	  	if (cell->prox.z < -NOISE) freebro = 1;
		cell->prox.x = cell->prox.y = cell->prox.z = 999.9;
		if (sis && cell->level == sis->level)       cell->prox.x = vdif(cell->xyz,cell->sis->xyz);
		if (sis && bro && sis->level == bro->level) cell->prox.y = vdif(cell->sis->xyz,cell->bro->xyz);
		if (bro && cell->level == bro->level)       cell->prox.z = vdif(cell->xyz,cell->bro->xyz);
		if (level == Data::depth) { // force ideal bond length for atom level
			cell->prox.x = bond;
			if (cell->prox.y>999.0) cell->prox.y = bond*1.85; // for 135 deg. angle
			cell->prox.z = bond;
		}
		if (freesis) { // mark sister bond as free
			cell->prox.x *= -1.0;
		}
		if (freebro) { // mark brother bond as free
			cell->prox.z *= -1.0;
		}
	}
	//Pi(level) Pi(cell->uid) Pi(cell->type) Pi(cell->sort) Pv(cell->geom) Pv(cell->dist) Pv(cell->prox) NL
	if (n == 0) return;
	for (i=0; i<n; i++) setAngle(cell->child[i], update);
}

void readlinks ( char *line, int maxin ) 
{
Cell	*ci, *cj;
FILE	*lin;
char	mode[5];
int	i,j, in, nlinks, resids = 0, at = 12;
int	linksin, print = 1;
//	read line: [LINKS|BONDS] [local|global] <filename> (local = internal numbering, global = PDB resid numbering)
	if (toupper(line[0])=='L') strcpy(mode,"Link"); else strcpy(mode,"Bond");
	if (islower(line[0])) print = 0;
	resids = 0;
	if (toupper(line[6])=='P') {	// pdbid
		resids = 1;
		if (print) printf(" # %sing by PDB resid numbers\n",mode);
	}
	if (toupper(line[6])=='L') {	// local
		resids = 2;
		if (print) printf(" # %sing by local sequential numbering\n",mode);
	}
	if (toupper(line[6])=='G') {	// global (uid)
		resids = 3; at++;
		if (print) printf(" # %sing by global sequential numbering\n",mode);
	}
	if (resids==0) { Pt(No valid numbering scheme specified\n)  return; }
	if (print) printf(" # reading: %s\n", line+at);
	lin = fopen(line+at,"r");
	linksin = 0;
	while (1) { int type, link, err = 0;
		in = read_line(lin,line);
		if (line[0]=='#') { Pt(#) Ps(line+1) NL continue; }
		if (in < 3) { fclose(lin); return; }
		sscanf(line,"%d %d %d %d", &i, &j, &type, &link); // type = thickness, link = end-end
		if (resids==1) { // PDB id (apply to all sequential pairs with those id.s)
			ci = cj = 0;
			DO1(n,Cell::total) { Cell *cn = Cell::uid2cell[n];
				if (cn->resn == i) ci = cn;
				if (cn->resn == j) cj = cn;
				if (ci && cj) { // got both so link, reset and try again
					if (mode[0]=='L') err = makeLink(ci,cj,type,link,print);
					if (mode[0]=='B') err = makeBond(ci,cj,type,link,print);
					if (err) { Ps(line) Pi(resids) NL exit(1); }
					ci = cj = 0;
				}
			}
		} else {
			if (resids==2) { // local
				i = atom2atom[i]; j = atom2atom[j];
				ci = atom2cell[i];
				cj = atom2cell[j];
			}
			if (resids==3) { // global
				ci = Cell::uid2cell[i];
				cj = Cell::uid2cell[j];
			}
			if (mode[0]=='L') err = makeLink(ci,cj,type,link,print);
			if (mode[0]=='B') err = makeBond(ci,cj,type,link,print);
			if (err) { Pt(#) Ps(line) Pi(resids) NL exit(1); }
		}
		linksin++;
		if (linksin == maxin) break;
	}
	fclose(lin);
}

int makeLink ( Cell *ci, Cell *cj, int type, int link, int print ) {
int	nlinks, m = ci->model;
	links = Data::model[m].links;
	nlinks = links[ci->level];
	if (nlinks>0) { int made = 0;
		for (int k=0; k<nlinks; k++) {
			if (ci->link[k].to == cj) {
				if (print) { Pt(# Existing link to) Pi(cj->uid) NL }
				return 0;
			}
		}
		for (int k=0; k<nlinks; k++) {
			if (ci->link[k].to == 0) { // free
				ci->link[k].to = cj;
				ci->link[k].type = type;
				ci->link[k].link = link;
				ci->link[k].next = link;
				// next is not needed in linking so use it to store extension limit 
				// 	(as % of length) with -ve link = unbreakable
				made = 1;
				break;
			}
		}
		if (made) {
			printf(" # Cell %d linked to cell %d   PDBids = %d to %d   dist = %f\n",
				ci->uid, cj->uid, ci->resn, cj->resn, vdif(ci->xyz,cj->xyz));
		} else {
			printf("NB: cell has no free links\n");
			return 1;
		}
	} else {
		printf("NB: tried to link a cell with no links\n");
		return 1;
	}
	return 0;
}

int makeBond ( Cell *ci, Cell *cj, int type, int link, int print ) {
int	nbonds, bi,bj, m;
	m = ci->model;
	bonds = Data::model[m].bonds;
	chain = Data::model[m].chain;
	nbonds = abs(chain[ci->level]);
	bi = -1;
	if (nbonds>0) {
		for (int k=0; k<=nbonds; k++) { // NB <=
			if (ci->bond[k].to == cj) {
				if (print) { Pt(# Existing bond to) Pi(cj->uid) NL }
				return 0;
			}
		}
		for (int k=0; k<=nbonds; k++) { // NB <=
			if (ci->bond[k].to == 0) { // free
				bi = k; break;
			}
		}
		if (bi < 0) {
			printf("NB: cell has no free bonds (uid=%d)\n",ci->uid); Pi(nbonds) NL
			return 1;
		}
	} else {
		printf("NB: tried to bond a cell with no bonds (uid=%d)\n",ci->uid);
		return 1;
	}
	m = cj->model;
	bonds = Data::model[m].bonds;
	nbonds = abs(chain[cj->level]);
	bj = -1;
	if (nbonds>0) {
		for (int k=0; k<=nbonds; k++) { // NB <=
			if (cj->bond[k].to == 0) { // free
				bj = k; break;
			}
		}
		if (bj < 0) {
			printf("NB: cell has no free bonds (uid=%d)\n",cj->uid); Pi(nbonds) NL
			return 1;
		}
	} else {
		printf("NB: tried to bond a cell with no bonds (uid=%d)\n",cj->uid);
		return 1;
	}
	ci->bond[bi].to   = cj;	cj->bond[bj].to   = ci;
	ci->bond[bi].next = bj;	cj->bond[bj].next = bi;
	ci->bond[bi].type =     cj->bond[bj].type = type;
	ci->bond[bi].link =     cj->bond[bj].link = link; // -ve flags a cross-link (not a chain-link)
	printf(" # Cell %d bonded to cell %d (link-type = %d)  PDBids = %d to %d   dist = %f\n",
		ci->uid, cj->uid, link, ci->resn, cj->resn, vdif(ci->xyz,cj->xyz));
	return 0;
}

void rebond () {
	Pi(nrelinks) NL
	FOR(i,nrelinks) // rewire from the atomic level
	{ Cell *cell = relink[i].at, *link = relink[i].to, *scope = relink[i].scope;
		if (relink[i].ends == 0) {	// terminus (cell = new N end, link = new C end)
			link->bond[1].to = link->bro = link->parent; link->bond[1].next = -1;
			cell->bond[0].to = cell->sis = cell->parent; cell->bond[0].next = -1;
			rewire(cell);
		}
		if (relink[i].ends == 1) {	// reconnect
			cell->bro = link;
			if (cell->bond) cell->bond[1].to = link;
			link->sis = cell;
			if (link->bond) link->bond[0].to = cell;
		}
		if (relink[i].ends == 2) {	// flag for unrefined bond length
			cell->prox.z = -999.9;	// -ve z = link to bro not refined
			link->prox.x = -999.9;	// -ve x = link to sis not refined
			cell->bond[1].link = link->bond[0].link = -999; // dito
		}
	}
	FOR(i,nrelinks) // upgrade/make bonds into hinges
	{ Cell *ca = relink[i].at, *cb = relink[i].to;
	  int	seta, setb, lin = -relink[i].ends, level = ca->level,
		nca = ca->nbonds+1, ncb = cb->nbonds+1; // always +1 for sis/bro
		if (lin <= 0) continue;
		seta = -1;
		FOR(i,nca) { if (ca->bond[i].to==cb) { seta = i; break; }}
		if (seta<0) { // not found so make new
			FOR(i,nca) { Cell *bond = ca->bond[i].to;
				if (bond==0 || bond->level<level) {
					seta = i; break;
				}
			}
			if (seta >= 0) ca->bond[seta].to = cb;
		}
		if (seta<0) {
			Pt(No free bond for HINGE A-side) NL
			Pi(ca->uid) Pi(ca->level) Pi(ca->type) NL
			exit(1);
		}
		setb = -1;
		FOR(i,ncb) { if (cb->bond[i].to==ca) { setb = i; break; }}
		if (setb<0) { // not found so make new
			FOR(i,ncb) { if (cb->bond[i].to==0) { setb = i; break; }}
			FOR(i,ncb) { Cell *bond = cb->bond[i].to;
				if (bond==0 || bond->level<level) {
					setb = i; break;
				}
			}
			if (setb >= 0) cb->bond[setb].to = ca;
		}
		if (setb<0) {
			Pt(No free bond for HINGE B-side) NL
			Pi(cb->uid) Pi(cb->level) Pi(cb->type) NL
			exit(1);
		}
		if (lin==1) { ca->bond[seta].link = 1; cb->bond[setb].link = 1; } // aN->bN, bN->aN
		if (lin==2) { ca->bond[seta].link = 2; cb->bond[setb].link = 3; } // aN->bC, bC->aN
		if (lin==3) { ca->bond[seta].link = 3; cb->bond[setb].link = 2; } // aC->bN, bN->aC
		if (lin==4) { ca->bond[seta].link = 4; cb->bond[setb].link = 4; } // aC->bC, bC->aC
	}
}

int follow ( Cell*, int );

int follow ( Cell *here, int to ) {
Bonds	*link = here->bond+to;
Cell	*into = link->to;
int	at = link->next;
int	mode;
	if (at<0) { // switch to sister path
		Pi(here->uid) Pt(end) NL
		return -1;
	}
	Pi(here->uid) Pi(into->uid) Pi(at) NL
	mode = follow(into,at);
	if (mode<0) {
		Pi(here->uid) Pi(here->nbro) NL
		if (here->nbro < 2) return -1;
		follow(here->bond[2].to,1);
	}
}

void rewire ( Cell *start ) {
int	i, j, k, n;
int	path[MAXIN];
float	d;
int	m = start->model,
	moltype = Data::model[m].moltype,
	subtype = Data::model[m].subtype,
	dna = moltype*subtype;
/*
Cell *c = start;
Cell::world->print();
LOOP{
Pi(c->uid) Pi(c->resn) Pi(c->parent->uid) NL
c = c->bro;
if (c->level<4) break;
}
Pi(start->level) Pi(start->uid) Pi(start->resn) NL
FOR(i,Cell::total) { Cell *c = Cell::uid2cell[i];
if (c->level>1 && c->uid==c->sis->sis->uid) { Pi(c->level) Pi(c->uid) Pi(c->sis->uid) Pi(c->sis->sis->uid) NL exit(1); }
}
{ Cell *c = Cell::uid2cell[18]; Pi(c->resn) Pi(c->sis->resn) Pi(c->sis->uid) Pi(c->bro->uid) Pi(c->bro->level) NL}
{ Cell *c = Cell::uid2cell[7];
Pi(c->kids) NL FOR(i,c->kids){ Cell *ci = c->child[i];
Pi(ci) NL
Pi(ci->resn) Pi(ci->sis->resn) Pi(ci->sis->uid) Pi(ci->bro->uid) Pi(ci->bro->level) NL}
NL }
Pi(Cell::uid2cell[7]->child[1]) NL
*/
	// propagate rewiring at the atomic level up through the higher levels 
	if (Data::depth < 2) return; // nothing to check (the world cannot be in a chain)
	for (i=Data::depth; i>1; i--)
	{ Cell	*lo = start, *hi = start->parent, *newlo, *newhi, *oldlo, *oldhi;
	  int	lolev = i, hilev = i-1,
		locha = Data::model[lo->model].chain[lolev],
		hicha = Data::model[hi->model].chain[hilev];
		if (locha < 1) continue;	// no chain at this level
		if (hicha < 1) break;	// no chain at this level
		for (j=0; j<hi->parent->kids; j++) { Cell *c = hi->parent->child[j];
			m = c->model;
			if (Data::model[m].moltype != moltype) continue; // don't mix models
			for (k=0; k<hicha; k++) { c->bond[k].to = 0; c->bond[k].next = -1; }
		}
		n = 0; path[0] = hi->uid;
		DO {
			if (lo->done) {
				Pt(Loop in new bond path:) NL
				Pi(lo->uid) Pi(lo->atom) Pi(lo->btom) Pi(lo->ctom) Pi(lo->level) NL
				exit(1);
			}
			lo->done = 1;
			if (path[n] != hi->uid) path[++n] = hi->uid;
			oldlo = lo->sis; oldhi = hi->sis;	// short for previous chain positions
			newlo = lo->bro; newhi = newlo->parent;	// next cell in low-level chain
			if (newlo->level != lolev) {	// end of the low-level chain
				if (hicha == 0) {
					Pt(# NB end of low chain in detached parent) NL
					continue; // no chain at high level
				} else { // end of lo chain so newlo is at hi level
					Pt(# NB end of low chain in chained parent) NL
					//if (moltype==1 && lolev==depth) {
					if (dna==1 && lolev==depth) { // ???????????
						// a nucleic acid chain ends at the start so don't mark as an end
						newlo->bro = newlo->bond[1].to = newlo; // reset to old brother
						n++;
					} else {
						newlo->bro = newlo->bond[1].to = newhi; // set end of parent chain
					}
					//newlo->bond[1].to = start->parent; // cyclic bond path
					//newlo->bond[1].next = 1;
					break;	  // must be end of level
				}
			}
			if (newhi != hi) {	// change of parent
				if (newhi == hi->bro) { // normal move into next high-level cell
					//Pt(normal move in parental chain) NL
				} else {	// the new parent is unexpected 
					//Pt(unexpected jump in parental chain) NL
					if (newhi == hi->sis) { // just been there
						//Pt(back to last) NL
					} else {
						//Pt(move to new parent) NL
					}
					hi->bro = newhi;
					newhi->sis = hi;
				}
			}
			lo = newlo; hi = newhi;	// next cell and parent
		}
		if (n==0) return;
		// path[] has collected the sequence of hi-level cells visited, now set bonds
		// bond[0]=back (->sis path), bond[1+]=out (->bro path)
		for (j=0; j<n; j++) { Cell *cj = Cell::uid2cell[path[j]]; // wipe old pointers
			cj->sis = cj->bro = Cell::world;
			cj->bond[0].to = cj->bond[1].to = Cell::world;
			cj->bond[0].next = cj->bond[1].next = -1;
			cj->nsis = cj->nbro = 0; // counter for bro links
		}
		for (j=0; j<n; j++)
		{ int	p = path[j], q = path[j+1], end = 0;
		  Cell	*last, *here = Cell::uid2cell[p], *next = Cell::uid2cell[q];
		  Bonds	*sis = here->bond+0, *bro = here->bond+(here->nbro+1);
			if (j==0) last = here->parent;
			if (next->bro==here) { // returning
				if (next==last) { // dead-end
					bro->to = here->bro = here->parent;
					bro->next = -1; // flag dead-end
					here->nbro++;
					sis->to = here->sis = last;
					sis->next = 0; // follow to next sis
				}
			} else { // set bro to next cell
				bro->to = here->bro = next;
				bro->next = 1; // follow to next bro
				here->nbro++;
				if (sis->next < 0) { // unset
					sis->to = here->sis = last;
					sis->next = 0; // follow to next sis
				}
			}
			last = here;
		}
/*
{ Cell *here = start = Cell::uid2cell[7]; int it;
DO {
Cell  *half = Cell::world;
Bonds *link, *sis = here->bond, *bro = here->bond+1;
Pi(it) Pi(here->uid) Pi(sis->to->uid) Pi(bro->to->uid) NL
if (it==0 && here==start) exit(1);
	it = 0;
	DO1(b,here->nbro) { // check available bro paths
		link = here->bond+b;
		if (link->to->nsis==0) { // new path
			if (link->next < 0) { // deadend
				here = sis->to;
				it = -1;
			} else { // found an open bro path
				here = link->to;
				here->nsis = 1;
				it = 1;
			}
			break;
		} 
		if (link->to->sort==1 && link->to->nsis==1) { // half open ladder
			half = link->to;
		}
	}
	if (it<0) { // deadend must retrace
		it = 0;
		continue;
	}
	if (it==1) continue;
	if (it==0) { // no open bro so follow sis, unless...
		if (here->sort==0 && bro->to->sort==1) { // move from loop into a stem
			here = bro->to;
			here->nsis = 2;
			continue;
		}
		if (half->uid > 0) {
			here = half;
			here->nsis = 2;
			continue;
		}
		if (here->sort==1) here->nsis = 2;
		here = sis->to;
		here->nsis = 2;
	}
}}
exit(1);
*/
		if (dna && i==depth ) // use half double stranded nucleic chain
		{ Cell	*last = Cell::uid2cell[path[n-1]],
			*past = Cell::uid2cell[path[n-2]],
			*half = Cell::uid2cell[path[n/2]];
				half->bro = half->parent;
				last->bro = past;
		} else {
			Cell::uid2cell[path[n-1]]->bro = Cell::uid2cell[path[n-1]]->parent;
		}
		start = start->parent;
		//
		/* print the path 
		Pt(Rewired path at level) Pi(i) NL
		{ Bonds *b = start->bond+1; int first = b->next;
		  Cell	*last = start;
			Pi(start->uid) Pi(first) NL
			LOOP{ // if b==last->bond the link is on the return path (sis) used for drawing bonds
				if (b==last->bond) { Pi(last->uid) Pt(-->) Pi(b->to->uid) Pi(last->bond[0].link) }
				if (b!=last->bond) Pi(b->to->uid) Pi(b->to->sis->uid) Pi(b->to->bro->uid) NL
				if (b->to->sis == b->to->bro) break; // dead-end
				last = b->to;
				b = b->to->bond+b->next;
				if (b->to==start && b->next==first) {
					Pi(last->uid) Pt(-->) Pi(b->to->uid) Pi(last->bond[0].link) NL
					break;
				}
			}
		}
		*/
	}
}

void setHinge ( char *line )
{
// HINGE <level> <id1> <id2> <type>
// or using UID:
// HINGE global <id1> <id2> <type>
//
int	n, a, b, lev, len, lin;		// <id[12]> = sequential number in <level>
Cell	*ca, *cb;			// <lin> = 1..4 NN,NC,CN,CC
int	seta, setb;
	if (line[6]=='g') { // use uid
		sscanf(line+12,"%d %d %d", &a, &b, &lin);
		ca = Cell::uid2cell[a];
		cb = Cell::uid2cell[b];
		lev = ca->level;
		if (lev != cb->level) {
			Pt(HINGE level mismatch) Pi(ca->level) Pi(cb->level) NL
			exit(1);
		}
	} else { // find id at given level
		sscanf(line+5,"%d %d %d %d", &lev, &a, &b, &lin);
		ca = cb = 0;
		n = 0;
		FOR(i,Cell::total) { Cell *cc = Cell::uid2cell[i];
			if (cc->level != lev) continue;
			if (n == a) ca = cc;
			if (n == b) cb = cc;
			n++;
		}
	}
	if (ca==0) { Pi(a) Pt(not found for HINGE) NL exit(1); }
	if (cb==0) { Pi(b) Pt(not found for HINGE) NL exit(1); }
	printf(" # HINGE set between level %d objects %d and %d of type %d\n", lev, a,b, lin);
	printf(" # HINGE cells are %d and %d types %d %d\n", ca->uid,cb->uid, ca->type,cb->type);
	if (ca->nbonds==0 || cb->nbonds==0) {
		Pt(No bonds allocated to make the HINGE) NL 
		Pi(ca->uid) Pi(ca->nbonds) NL
		Pi(cb->uid) Pi(cb->nbonds) NL
		exit(1);
	}
	relink[nrelinks].at = ca;
	relink[nrelinks].to = cb;
	relink[nrelinks].ends = -lin; // -ve flag = convertion-to/addition-of a hinge
	nrelinks++;
}
