#include "main/common.h"

pthread_mutex_t linker_lock;
pthread_mutex_t drawBonds_lock;
int argcin; char** argvin;

void driver();
void bumper();
void shaker();
void keeper();
void viewer();

float *focus;
float scalein, scaleout;

//           100000000 = 1/10 sec.
long nanos = 100000000;	// good speed for watching
int  delay = 10;//10;

float	atom2pdb;	// scale factor

char	line[222];
int	teratoms, endatoms, allatoms, atoms, secs;
Vec	cent;

Cells	**atom1cell, **atom2cell, **secs2cell;
int	 *atom2atom,  *atom2resn,  *resn2atom;

typedef struct {
	Cells	*at, *to, *scope; // cell source, target, scope of relink
	int	ends;	  // 0 = relink, -1 = new start, 1 = new finish
} Relinks;
Relinks	*relink;
int	nrelinks;

typedef struct {
	Cells	*a, *b;
	float	d;
} Pairs;

sortAlpha ( const void *ac, const void *bc ) {
Pairs   *a = (Pairs*)ac, *b = (Pairs*)bc;
	if (a->d < b->d) return -1;
	if (a->d > b->d) return  1;
	return 0;
}

void findPole ( Vec, Vec, Vec, float, float, float, Vec* );
void reset2axis ( Cells*, int, float );

int i;

void *acts ( void *ptr )
{
struct  timespec wait;
char	*message;
	message = (char *) ptr;
	sleep(1); // delay start to allow ends to be set
	printf("Starting %s \n", message);
	wait.tv_sec = 0;
	wait.tv_nsec = nanos;
	while (1) {
		if (data[0] > 0 && data[5] == 999 ) {
			// give time for presort
			printf("Waiting for presort\n");
			if (data[1] < 10000) sleep(1); else sleep(10);
			continue;
		}
		nanosleep(&wait,NULL);
		if (delay > 0) driver(data,world);
		data[0]++;
		if (*message=='D') break;
	}
	printf("Stopping %s \n", message);
}

void *fixs ( void *ptr )
{
float	weight = (float)nanos/100000000.0; // = 1 when speed = 1, .1 for 10, .01 for 100
struct  timespec wait;
char	*message;
	message = (char *) ptr;
	sleep(1); // delay start to allow ends to be set
	printf("Starting %s \n", message);
	weight *= 0.1;
	wait.tv_sec = 0;
	wait.tv_nsec = nanos;
	while (1) {
		fixSSEaxis(world,0,weight);
		nanosleep(&wait,NULL);
/*
		fixRNAaxis(world,0,weight);
		nanosleep(&wait,NULL);
		fixRNAstem(world,0,weight);
		nanosleep(&wait,NULL);
		fixDOMaxis(world,0,1.0);
*/
		fixANYaxis(world,0,weight);
		nanosleep(&wait,NULL);
		reset2axis(world,0,weight);
		if (*message=='D') break;
	}
	printf("Stopping %s \n", message);
}

void *move ( void *ptr )
{
struct  timespec wait;
char	*message;
	message = (char *) ptr;
	printf("Starting %s \n", message);
	wait.tv_sec = 0;
	wait.tv_nsec = nanos;
	while (1) {
		shaker(data,world);
		nanosleep(&wait,NULL);
//		keeper(data,world,types);
		if (*message=='D') break;
	}
	printf("Stopping %s \n", message);
}

void *link ( void *ptr )
{
struct  timespec wait;
char	*message;
	message = (char *) ptr;
	printf("Starting %s \n", message);
	wait.tv_sec = 0;
	wait.tv_nsec = nanos;
	while (1) {
		nanosleep(&wait,NULL);
		linker(data,world,types);
		if (*message=='D') break;
	}
	printf("Stopping %s \n", message);
}

void *rank ( void *ptr )
{
long	sleep;
struct  timespec wait;
int	bumping;
char	*message;
	message = (char*)ptr;
	printf("Starting %s \n", message);
	wait.tv_sec = 1; // <------------NB + 1 sec.
	wait.tv_nsec = nanos;
	while (1) {
		sorter(data,world,scene);
		nanosleep(&wait,NULL);
		if (*message=='D') break;
	}
	printf("Stopping %s \n", message);
}

void *bump ( void *ptr )
{
long	sleep;
struct  timespec wait;
int	bumping;
char	*message;
	message = (char*)ptr;
	printf("Starting %s \n", message);
	sscanf(message+7,"%d", &bumping);
	wait.tv_sec = 0;
	wait.tv_nsec = nanos;
	while (1) {
		bumper(data,world,scene,bumping);
		nanosleep(&wait,NULL);
		if (*message=='D') break;
	}
	printf("Stopping %s \n", message);
}

void *view ( void *ptr )
{
	char *message;
	message = (char *) ptr;
	printf("Starting %s \n", message);
	viewer(focus,data,world,scene,types,argcin,argvin);
	printf("Stopping %s \n", message);
}

void *rest ( void *ptr )
{
struct  timespec wait;
char	*message;
	message = (char *) ptr;
	printf("Starting %s \n", message);
	wait.tv_sec = 0;
	wait.tv_nsec = nanos;
	sorter(data,world,scene);
	data[0] = 0;
	while (1) {
		Pt(bumper) NL
		bumper(data,world,scene,0);
		Pt(driver) NL
		driver(data,world);
		Pt(shaker) NL
		shaker(data,world);
		Pt(keeper) NL
		keeper(data,world,types);
		Pt(linker) NL
		linker(data,world,types);
		Pt(sorter) NL
		sorter(data,world,scene);
		data[0]++;
		nanosleep(&wait,NULL);
	}
	printf("Stopping %s \n", message);
}

//define MAXIN 100000 (needed for 3chy/unit3+.run)
#define  MAXIN 10000

void startup ( FILE *pdb )
{
int	i,j,k,n, wait;
char	line[222];
	Pi(nmodels) NL
	Escale = ESCALE; // ellipsoid range factor
	// allocate pairs of cell to collect chain relinks in pass1set
	relink = (Relinks*)alloca(sizeof(Relinks)*9999);
	nrelinks = 0;
	// set-up the type structure
	types = (Types***)malloc(sizeof(Types**)*nmodels); TEST(types)
	for (i=0; i<nmodels; i++) types[i] = 0;
	// set-up the cell structure
	world = (Cells*)malloc(sizeof(Cells)); TEST(world)
	secs = teratoms = endatoms = allatoms = total = 0;
	uid2cell  = (Cells**)malloc(sizeof(Cells*)*MAXIN*8); TEST(uid2cell) // global back-ref.lookup (total)
	atom2cell = (Cells**)alloca(sizeof(Cells*)*MAXIN*6);	// lookup array for all atoms (allatoms)
	secs2cell = (Cells**)alloca(sizeof(Cells*)*MAXIN*1);	// WAS 10000);	// lookup for all SSEs
	atom1cell = (Cells**)alloca(sizeof(Cells*)*MAXIN*1);	// WAS 1000);	// lookup for domain (teratoms)
	atom2resn = (int*)alloca(sizeof(int)*MAXIN*6);	// lookup from global atom id to pdb resid (allatoms)
	atom2atom = (int*)alloca(sizeof(int)*MAXIN*4);	// lookup from local to global atom id (endatoms)
	resn2atom = (int*)alloca(sizeof(int)*MAXIN*1);	// lookup from global pdb resid to atom id
	world->type = 0;
	world->model = 0;
	vinit(&(world->xyz));
	vinit(&(world->endN));
	vinit(&(world->endC));
	vinit(&(world->cent));
	world->endN.z =  1.0;
	world->endC.z = -1.0;
	vinit(&cent);
	pass1set(world,0,-1,pdb,0);
	vdiv(&cent, (float)allatoms);
	Pi(allatoms) Pv(cent) NL
	pass2set(world,0);
	pass3set(world,0);
	if (nrelinks) rebond();
	// read cross-link data
	readlinks(line,pdb); // check to end of file for LINKS
	fclose(pdb);
	// set-up the scene data
	world->solid = -1;
	world->lost = 0;
	world->bump = 0;
	world->done = 0;
	world->live = 999;
	for (i=0; i<SPARE; i++) { // contents allocated in copyCell()
		n = world->kids + i;
		world->child[n]->kids = -1; // flag that no childern exist
		world->child[n]->nbonds = -1; // or bonds to them
		world->child[n]->nlinks = -1; // or links to them
	}
	world->far = scaleout;	// world is never moved so using far to hold output scale factor
	world->link = uid2cell; // world is never linked so using link to pass uid2cell address
	scene = (Cells*)malloc(sizeof(Cells)); TEST(scene)
	scene->rank = (Cells***)malloc(sizeof(Cells**)*3); TEST(scene->rank)
	for (j=0; j<3; j++) {
		scene->rank[j] = (Cells**)malloc(sizeof(Cells*)*total*2); TEST(scene->rank[j])
	}
	total = 0;
	sceneCell(world,0);	// sets-up front-view
	sceneCell(world,0);	// sets-up rearview
	scene->far = shrink;	// shrink factor applied to ever cell->xyz on every cycle
	total /= 2;
	wait = (int)cbrt((float)total)*2;
	data[0] = -wait;	// frame counter (-ve = settle time before presort)
	data[1] = total;	// total bodies
	data[2] = depth;	// deepest level
	data[3] = 5000;		// distance cooling factor
	data[4] = nmodels;	// number of parameter sets 
	data[5] = 999;		// view sort 123=XYZ, +/-=front/back (999 = flag to pre-sort lists)
	data[6] = 10;		// freezing point
	data[7] = 0;		// used for default moltype (multiple moltypes no implemented yet)
	data[8] = allatoms;	// number of atoms (lowest level)
	data[9] = -1;		// flag for types set in viewer() -1=unset, then copyCell() update flag
	Pi(total) Pi(depth) Pi(allatoms) Pi(secs) NL
	printf("Data structures all setup\n");
}

pass1set ( Cells *cell, int level, int id, FILE *file, int subfile )
// first pass allocates memory, sets values and gets atoms
{
FILE	*pdb;
float	x,y,z, s;
int	i,j,k,m,n, in, sort = 0;
int	nlinks, nbonds, newfile = 0;
Vec	heli, axis, grow, tran, xrot,yrot,zrot;
float	helix, theta, spinx, spiny, spinz;
char	nextfile[111];
	helix = theta = spinx = spiny = spinz = 0.0;
	vinit(&axis); vinit(&grow); vinit(&tran);
	vinit(&xrot); vinit(&yrot); vinit(&zrot);
	cell->atom = cell->ends = 0;
	in = read_line(file,line);
	if (level) vcopy(cell->parent->xyz, &(cell->xyz));
	while (line[0]=='#') {	// comments
		Ps(line) NL
		in = read_line(file,line);
	}
	if (line[0]=='M') { // read new MODEL type
		NL Ps(line) NL
		sscanf(line+6,"%d", &model);
		printf("\nSwitching to model %d\n", model);
		m = model*M*N+N; // M rows of N columns for each parameter set
		// first row of data is kept for global parameters
		moltype = data[m]; // first position specifies mol.type: 0=prot, 1=RNA
		if (moltype==0) printf("Molecule type = PROTein\n");
		if (moltype==1) printf("Molecule type = RNA\n");
		if (moltype==2) printf("Molecule type = CHEMical\n");
		names = data+m+N*1; class = data+m+N*2;
        	sizes = data+m+N*3; bumps = data+m+N*4; links = data+m+N*5; chain = data+m+N*6;
		kicks = data+m+N*7; keeps = data+m+N*8; repel = data+m+N*9; bonds = data+m+N*10;
		split = data+m+N*11; local = data+m+N*12;
		in = read_line(file,line);
	}
	s = scalein;
	nbonds = abs(chain[level]);
	nlinks = links[level];
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
		vcopy(cell->xyz, &axis);
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
	vinit(&(cell->endN)); vinit(&(cell->endC));
	cell->endN.z = cell->endC.z = 1234.5;	// marker for unset ends
	cell->cent.x = cell->cent.y = cell->cent.z = -1.0; // -ve dist = unset
	if (line[0] == 'G') {	// GROUP <sort> <members> [<x> <y> <z>  <x> <y> <z> ]
		sscanf(line+6,"%d %d", &sort, &n);
		if (in > 20) { float x1,y1,z1, x2,y2,z2;
			sscanf(line+11,"%f %f %f  %f %f %f", &x1,&y1,&z1, &x2,&y2,&z2);
			cell->endN.x = x1*s; cell->endN.y = y1*s; cell->endN.z = z1*s;
			vsum(cell->xyz, &(cell->endN));
			cell->endC.x = x2*s; cell->endC.y = y2*s; cell->endC.z = z2*s;
			vsum(cell->xyz, &(cell->endC));
			cell->ends = 1;
		}
		if (level == depth-1) { // SSE level
			secs2cell[secs] = cell;
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
		vsum(cell->xyz, &cent);
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
		atom1cell[teratoms] = cell;	// local to current file
		atom2cell[allatoms] = cell;	// global over all atoms
		atom2atom[endatoms] = allatoms; // relink scope to global
		sscanf(line+22,"%d", &resn);
		cell->resn = resn;
		if (resn<0) { Pt(*NB* negative residue numbers are not good: default to 0) NL resn = 0; }
		resn2atom[resn] = allatoms;
		atom2resn[allatoms] = resn;
	}
	cell->uid = total;
	uid2cell[total] = cell;
	total++;
	if (total > MAXIN*8) { Pt(Increase MAXIN) Pi(total) NL exit(1); }
	if (n==1 && chain[level+1]) { // an only child is not a chain (but might be later)
		printf("*NB* chain of one at level %d\n", level+1);
		// chain[level+1] = 0; 
	}
	cell->model = model;
	cell->kids = cell->cots = n;
	cell->id = id;
	cell->solid = 0;
	cell->far = 9.9;
	for (i=0; i<NDATA; i++) cell->data[i] = 0;
	cell->live = 9; cell->done = cell->lost = cell->bump = cell->busy = cell->empty = 0;
	cell->level = level;
	cell->type = class[level];	// pick-up the type assigned to the current level
	//if (cell->type==3) sort = E/2;  // default = 1:1 (sphere)
	cell->sort = sort;
	identM(cell->rot);
	cell->nbonds = nbonds;
	if (nbonds > 0) {
		if (nbonds > 9) { Pt(*NB*) Pi(nbonds) NL exit(1); }
		cell->bond = (Bonds*)malloc(sizeof(Bonds)*nbonds); TEST(cell->bond)
		for (i=0; i<nbonds; i++) {
			cell->bond[i].to = 0;	// bro=0, sis=1
			cell->bond[i].type = 0; // default is unset (0)
			cell->bond[i].next = 0; // default is follow brother
		}
	} else { cell->bond = 0; }
	cell->nlinks = nlinks;
	if (nlinks > 0) {
		cell->link = (Cells**)malloc(sizeof(Cells*)*nlinks); TEST(cell->link)
		for (i=0; i<nlinks; i++) cell->link[i] = 0;
	} else { cell->link = 0; }
	for (i=0; i<3; i++) cell->ranks[i] = id;
	cell->junior = cell->senior = 0;
	cell->starts = cell->finish = 0;
	cell->nsis = cell->nbro = 0;
	cell->sis = cell->bro = 0;
	if (n==0) return;
	// make ranked lists for the children
	if (level>0) m=n; else m=n+SPARE;  // make spare cells at world's end
	cell->rank = (Cells***)malloc(sizeof(Cells**)*m); TEST(cell->rank)
	for (i=0; i<m; i++) {
		cell->rank[i] = (Cells**)malloc(sizeof(Cells*)*3); TEST(cell->rank[i])
	}
	// make the children
	cell->child = (Cells**)malloc(sizeof(Cells*)*m); TEST(cell->child)
	for (i=0; i<n; i++) {
		cell->child[i] = (Cells*)malloc(sizeof(Cells)); TEST(cell->child[i])
		cell->child[i]->parent = cell;
		pass1set(cell->child[i],level+1,i,pdb,subfile);
		for (j=0; j<3; j++) cell->rank[i][j] = cell->child[i];
	}
	if (level==0) { // make but don't fill extra world child
		for (m=n; m<n+SPARE; m++) {
			cell->child[m] = (Cells*)malloc(sizeof(Cells)); TEST(cell->child[m])
			cell->child[m]->parent = world;
		}
	}
	// twists and turns implenented on return path (NB correct input and execution order is critical)
	// HELIX on Z axis: shift by x,y, spin theta on Z then shift z along Z
	if (helix*helix>NOISE && vsqr(heli)>NOISE) { float z = heli.z;
		printf("Cell %d at level %d helical rotation by %f rad.\n", cell->id, cell->level, helix);
		moveCell(cell,heli,1);
		vcopy(cell->xyz, &heli); heli.z += 1.0;
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
		vcopy(cell->xyz, &axis); axis.x += 1.0;
		spinCell(cell,cell->xyz,axis,spinx);
	}
	if (spiny*spiny > NOISE) {
		printf("Cell %d at level %d rotated around Y by %f rad.\n", cell->id, cell->level, spiny);
		vcopy(cell->xyz, &axis); axis.y += 1.0;
		spinCell(cell,cell->xyz,axis,spiny);
	}
	if (spinz*spinz > NOISE) {
		printf("Cell %d at level %d rotated around Z by %f rad.\n", cell->id, cell->level, spinz);
		vcopy(cell->xyz, &axis); axis.z += 1.0;
		spinCell(cell,cell->xyz,axis,spinz);
	}
	// TRANSlate (after spins)
	if (vsqr(tran) > NOISE) { float d = vmod(tran);
		printf("Cell %d at level %d moved by %f\n", cell->id, cell->level, d);
		moveCell(cell,tran,1);
	}
	// SCALE 
	if (vsqr(grow) > NOISE) growCell(cell,grow);
	if (newfile || level==0) {	// check for information at the end of the old file before new is opened
		in = read_line(pdb,line);
		if (teratoms) {
			if (line[0]=='S') {		// SHEET follows
				printf("Beta sheet links set\n");
				sheetset(pdb);		// read (to EOF or END) and set BETA pairs
			}
			in = read_line(pdb,line);
		}
		while (line[0]=='R')			// REBOND+RELINK
		{ int	at, to; Cells *cat, *cto;
			Ps(line) NL
			in = 12; if (toupper(line[7])=='G') in++;
			sscanf(line+in,"%d %d", &at, &to);
			if (toupper(line[7])=='P') { // use pdbid
				cat = atom2cell[resn2atom[at]];
				cto = atom2cell[resn2atom[to]];
			}
			if (toupper(line[7])=='L') { // local numbering (sequential from start, INPUT0 or REBOND)
				cat = atom2cell[atom2atom[at]]; // atom2atom[] gives the global (allatom) id 
				cto = atom2cell[atom2atom[to]];
			}
			if (toupper(line[7])=='G') { // global (uid) numbering
				cat = uid2cell[at];
				cto = uid2cell[to];
			}
			relink[nrelinks].at = cat;
			relink[nrelinks].to = cto;
			relink[nrelinks].ends = 1;	// flag to refine bond length
			relink[nrelinks].scope = cell;	// links external to the scope are preserved in relinks()
			if (line[2]=='L') relink[nrelinks].ends = 2; // flag to leave bond unrefined/unrendered
			nrelinks++;
			endatoms = 0;
			in = read_line(pdb,line);
			if (in<6) break;
		}
		while (line[0]=='E' && line[3]=='S')	// ENDS relocates the c-terminus
		{ int	at, to; Cells *cat, *cto;	// at=old, to=new, flag=0  NB fix needed for pdbid
			Ps(line) NL
			sscanf(line+5,"%d %d", &at, &to);
			cat = atom2cell[atom2atom[at]];
			cto = atom2cell[atom2atom[to]];
			relink[nrelinks].at = cat;
			relink[nrelinks].to = cto;
			relink[nrelinks].ends = 0;
			nrelinks++;
			endatoms = 0;
			teratoms = 0;
			in = read_line(pdb,line);
			if (in<6) break;
		}
		while (line[0]=='H' && line[1]=='I')	// HINGE <level> <id1> <id2> <type>
		{ int	n, a, b, lev, len, lin;		// <id[12]> = sequential number in <level>
		  Cells *ca, *cb, *cc;			// <len> = length*10, <lin> = 1..4 NN,NC,CN,CC
		  int	seta, setb, make = 5;
			Ps(line) NL
			sscanf(line+6,"%d %d %d %d %d", &lev, &a, &b, &len, &lin);
			ca = cb = 0;
			n = 0;
			for (i=1; i<total; i++) {
				cc = uid2cell[i];
				if (cc->level != lev) continue;
				if (n == a) ca = cc;
				if (n == b) cb = cc;
				n++;
			}
			if (ca==0) { Pi(a) Pt(not found for HINGE) NL exit(1); }
			if (cb==0) { Pi(b) Pt(not found for HINGE) NL exit(1); }
//Pi(ca->level) Pi(ca->id) Pi(ca->type) Pi(cb->level) Pi(cb->id) Pi(cb->type) NL
			// store hinge as CHEM bond (next = link, type = length)
			if (ca->bond==0) {
				ca->bond = (Bonds*)malloc(sizeof(Bonds)*make); TEST(ca->bond)
				for (i=0; i<make; i++) ca->bond[i].to = 0;
			}
			if (cb->bond==0) {
				cb->bond = (Bonds*)malloc(sizeof(Bonds)*make); TEST(cb->bond)
				for (i=0; i<make; i++) cb->bond[i].to = 0;
			}
			for (i=0; i<make; i++) { if (ca->bond[i].to==0) { seta = i; break; }}
			for (i=0; i<make; i++) { if (cb->bond[i].to==0) { setb = i; break; }}
			ca->bond[seta].to = cb; cb->bond[setb].to = ca; // two (a-b, b-a) links set
			if (lin==1) { ca->bond[seta].next = 1; cb->bond[setb].next = 1; } // aN->bN, bN->aN
			if (lin==2) { ca->bond[seta].next = 2; cb->bond[setb].next = 3; } // aN->bC, bC->aN
			if (lin==3) { ca->bond[seta].next = 3; cb->bond[setb].next = 2; } // aC->bN, bN->aC
			if (lin==4) { ca->bond[seta].next = 4; cb->bond[setb].next = 4; } // aC->bC, bC->aC
			ca->bond[seta].type = cb->bond[setb].type = len;
			in = read_line(pdb,line);
			if (in<6) break;
		}
		if (line[0]=='L' && line[4]=='S') {	// LINKS <mode> <filename>
			Ps(line) NL
			readlinks(line,0);
		}
		if (level) fclose(pdb);
		teratoms = 0;
		subfile--;
		return;
	}
}

pass2set ( Cells *cell, int level )
// second pass sets links and relationships
{
int	i,j,k, m,n, id, nlinks;
	id = cell->id;
	if (level > depth) depth = level;
	n = cell->kids;
	m = cell->model*M*N+N;
	moltype = data[m];
	links = data+m+N*5; chain = data+m+N*6;
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
	return;
}

getWcent ( Cells* cell, Vec* out )
// get the centroid of the children weighted by their number of kids
{
int	i;
Vec	cog, add;
float	wkids = 0.0;
	if (cell->kids<1) {
		out->x = cell->xyz.x;
		out->y = cell->xyz.y;
		out->z = cell->xyz.z;
		return;
	}
	vinit(&cog);
        for (i=0; i<cell->kids; i++)
        { float w = 1.0+(float)(cell->child[i]->kids); // weight into sum by kids
               	vcopy(cell->child[i]->xyz, &add);
               	vmul(&add,w); vsum(add, &cog);
               	wkids += w;
        }
        vdiv(&cog,wkids);            // find weighted centre of children
	out->x = cog.x;
	out->y = cog.y;
	out->z = cog.z;
}

setWcent ( Cells* cell ) {
Vec	cog;
	getWcent(cell, &cog); // not substituted directly to avoid vinit(cell->xyz)
	vcopy(cog, &(cell->xyz));
}

pass3set ( Cells *cell, int level )
// third pass adds links across chain of chains (can only be done after pass2set() )
// sums atoms into higher level centres and sets ends for ellipsoids and tubes
{
int	i,j,k, m,n, id, type, sort;
	id = cell->id;
	n = cell->kids;
	type = abs(cell->type);
	sort = cell->sort;
	m = cell->model*M*N+N;
	moltype = data[m];
	chain = data+m+N*6;
	sizes = data+m+N*3;
	if (n==0) {
		vsub(cell->xyz, cent, &(cell->xyz));
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
			cell->starts->bond[1].to = cell->starts->sis = cell->sis->finish;
		}
		if (cell->level == cell->bro->level) {	// set brother of last child as first child of brother
			cell->finish->bond[0].to = cell->finish->bro = cell->bro->starts;
		}
	}
	for (i=0; i<n; i++) pass3set(cell->child[i], level+1);
	if (level < depth) { // sum children into upper levels
		setWcent(cell);
		if (cell->ends) // ends have been read in (just to set the axis not its ends)
		{ float	size = 0.1*(float)sizes[level], // for tube and ellipsoids, size = section diameter
		  	axis = vdif(cell->endC,cell->endN),
			fn = (float)cell->kids;
		  Vec	x;
			vsub(cell->endC,cell->endN, &x);
			vnorm(&x); // unit axis in x
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
			// make symmetric about centre
			vmul(&x,axis*0.5);
			vsub(cell->xyz, x, &(cell->endN));
			vadd(cell->xyz, x, &(cell->endC));
			if (moltype==1 && type==2 && sort==1) // for RNA SSE put endC near termini
			{ float dnn, dcn, dnc, dcc;
				dnn = vdif(cell->endN,cell->junior->xyz);
				dnc = vdif(cell->endC,cell->junior->xyz);
				dcn = vdif(cell->endN,cell->senior->xyz);
				dcc = vdif(cell->endC,cell->senior->xyz);
				if (dnn+dcn < dnc+dcc) { // swap ends
					vsub(cell->xyz, x, &(cell->endC));
					vadd(cell->xyz, x, &(cell->endN));
				}
			}
			DO(i,n) { Cells *kidi = cell->child[i]; // set internal reference distances
				kidi->cent.x = vdif(kidi->xyz,cell->endN);
				kidi->cent.y = vdif(kidi->xyz,cell->xyz);
				kidi->cent.z = vdif(kidi->xyz,cell->endC);
			}	
		}
	}
}

rebond () {
int	i, j, k;
float	d;
	Pi(nrelinks) NL
	for (i=0; i<nrelinks; i++)
	{ Cells *cell = relink[i].at, *link = relink[i].to, *scope = relink[i].scope;
		if (relink[i].ends < 1) continue;
		if (relink[i].ends == 1) {
			cell->bro = link;
			cell->parent->bro = link->parent;
			link->sis = cell;
			link->parent->sis = cell->parent;
		}
		if (relink[i].ends == 2) {	// flag for unrefined bond length
			cell->prox.z = -999.9;	// -ve z = link to bro not refined
			link->prox.x = -999.9;	// -ve x = link to sis not refined
			cell->bro = cell->parent; // so won't be rendered
			link->sis = link->parent;
		}
	}
//exit(1);
}
/*
relinks () {
int	i, j, k;
	Pi(nrelinks) NL
	for (i=0; i<nrelinks; i++) // reset the bond counters (nsis, nbro)
	{ Cells *cell = relink[i].at, *link = relink[i].to, *scope = relink[i].scope;
		if (relink[i].ends < 1) continue;
		if (cell->bond==0) continue;
		while (1) { int	level = cell->level; Cells *elder;
			if (cell == link) break;
			if (cell->bro==0 || link->sis==0) break;
			// initialise forward link counters
			elder = cell->bro;
			for (j=level; j>scope->level; j--) {
				elder = elder->parent;
				if (elder == scope) { // cell's brother is a member of the scope family
					cell->nbro = 0;
					break;
				}
			}
			// initialise backward link counters
			elder = link->sis;
			for (j=level; j>scope->level; j--) {
				elder = elder->parent;
				if (elder == scope) { // link's sister is a member of the scope family
					link->nsis = 0;
					break;
				}
			}
			if (link->nsis && link->id<999) link->id += 999; // flag to keep this link in relinksis()
			cell = cell->parent; link = link->parent;
		}
	}
	for (i=0; i<nrelinks; i++) // process ENDS first before they get overwritten
	{ Cells *cell = relink[i].at, *newt = relink[i].to;
		if (relink[i].ends > 0) continue;
		if (cell->bond==0) continue;
		while (1) {
			cell->nbro = 0;
			newt->nbro = 1;
			newt->bond[1] = newt->bro = cell->bro; // set forward link
			newt->parent->finish = newt;
			cell = cell->parent;
			newt = newt->parent;
			if (cell->id < 0) break;
			if (cell->sis==0 || cell->bro==0) break;
		}
	}
	for (i=0; i<nrelinks; i++) // process relinks
	{ Cells *cell = relink[i].at, *link = relink[i].to;
		if (relink[i].ends < 1) continue;
		if (cell->bond==0) continue; // skip ends
		while (1) { int	forks = abs(chain[cell->level]);
			if (cell == link) break;
			if (cell->bro==0 || link->sis==0) break;
			if (relink[i].ends == 2) {	// flag for unrefined bond length
				cell->prox.z = -999.9;	// -ve z = link to bro not refined
				link->prox.x = -999.9;	// -ve x = link to sis not refined
			}
			if (cell->nbro < forks) {
				cell->bond[cell->nbro][1] = link;
				cell->nbro++;
			}
			if (link->nsis < forks) {
				link->bond[link->nsis][0] = cell;
				link->nsis++;
			}
			cell->bro = cell->bond[0][1];
			link->sis = link->bond[0][0];
			cell = cell->parent; link = link->parent;
		}
	}
	relinksis(world);
	printf("Relinked OK\n");
}

relinksis ( Cells *cell) {
int	i,j,k, m,n = cell->kids, level = cell->level;
	m = cell->model*M*N+N;
	chain = data+m+N*6;
	if (chain[level+1] > 1) { // children are branched so reconstruct the sisters
		for (i=0; i<n; i++) { Cells *ci = cell->child[i];	// for each cell
			ci->nsis = 0;					// find who links in
			if (ci->bond[0][0]->level==level) ci->nsis = 1; // keep upper level links
			if (ci->id >= 999) ci->nsis = 1;		// keep links out of scope
			for (j=0; j<n; j++) { Cells *cj = cell->child[j];
				for (k=0; k<cj->nbro; k++) { Cells *ck = cj->bond[k][1];
					if (ck != ci) continue;
					ci->bond[ci->nsis][0] = cj;
					ci->nsis++; break;
				}
			}
			ci->sis = ci->bond[0][0];
		}
	}
	for (i=0; i<n; i++) { Cells *ci = cell->child[i];
		relinksis(ci);
		if (ci->id >= 999) ci->id -= 999;
	}
}
*/

setangle ( Cells *cell, int level, int update )
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
int	i, n, id;
float	f, g, fix = 0.1;
Vec	alphA, alphD, betaA, betaD;
float	s = bondCA/3.8; // scale internal/real
	alphA.x = alphA.z = -0.9; alphA.y =  1.6;
	betaA.x = betaA.z =  2.7; betaA.y =  2.2;
	alphD.x = alphD.z = 5.1*s; alphD.y =  6.2*s;
	betaD.x = betaD.z = 9.9*s; betaD.y = 13.5*s;
//if(level==0) {Pr(bondCA) Pv(betaA) Pv(alphA) Pv(betaD) Pv(alphD) NL}
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
			if (cell->bro->bro) {
				if (angle(cell->xyz, cell->bro->xyz, cell->bro->bro->xyz) > NOISE) {
					tau2 = torsion(cell->bro->bro->xyz,cell->bro->xyz,cell->xyz,cell->sis->xyz);
				}
				dis2 = vdif(cell->bro->bro->xyz,cell->sis->xyz); // i+2..i-1
			}
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
	{ float bond = 0.1*(float)(sizes[level]+bonds[level]);
	  Cells	*sis = cell->sis, *bro = cell->bro;
	  int	freesis, freebro; // -ve dist = free bond (first set by RELINK in relinks())
	  	freesis = freebro = 0;
	  	if (cell->prox.x < -NOISE) freesis = 1;
	  	if (cell->prox.z < -NOISE) freebro = 1;
		cell->prox.x = cell->prox.y = cell->prox.z = 999.9;
		if (sis && cell->level == sis->level)       cell->prox.x = vdif(cell->xyz,cell->sis->xyz);
		if (sis && bro && sis->level == bro->level) cell->prox.y = vdif(cell->sis->xyz,cell->bro->xyz);
		if (bro && cell->level == bro->level)       cell->prox.z = vdif(cell->xyz,cell->bro->xyz);
		if (level == depth) { // force ideal bond length for atom level
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
//Pi(level) Pi(cell->id) Pi(cell->type) Pv(cell->geom) Pv(cell->dist) Pv(cell->prox) NL
	if (n == 0) return;
	for (i=0; i<n; i++) setangle(cell->child[i], level+1, update);
}

sheetset ( FILE *beta ) {
int	i;
	// set the beta sheet links
	while (1) { int io, a,b,c; Cells *from, *to;
		io = read_line(beta,line);
		if (io < 10) break;
		sscanf(line+6,"%d %d %d", &a, &b, &c);
		if (c>5) { Ps(line) NL }
		from = atom1cell[a]; to = atom1cell[b];	// set links at CA level (uses domain ids)
		//if (from->parent->sort<2 || to->parent->sort<2) continue; // SHEET atoms not in beta
		for (i=2; from->link && i<links[from->level]; i++) { // 0,1 used for +/-2 in strand
			if (from->link[i] == to) break; // already set
			if (from->link[i] == 0) {
				from->link[i] = to;
				break;
			}
		}
		from = from->parent; to = to->parent;	// set links at strand (midpoint) level
		if (from->sort<2 || to->sort<2) continue; // some SHEET2 atoms may be in loops
		for (i=0; from->link && i<links[from->level]; i++) {
			if (from->link[i] == to ) break;  // already set
			if (from->link[i] == 0) {
				from->link[i] = to;
				break;
			}
		}
	}
}

setlinks () {
int	i, j, m, n, hold = secs*10;
int	nlinks = 0;
Pairs	*pairs;
	// set-up the pairwise alpha-alpha, alpha-beta and loop-loop links
	bondCA = 0.1*(float)(sizes[depth]+bonds[depth]);
	/// pairs = (Pairs*)alloca(sizeof(Pairs)*hold); TEST(pairs) *NB* this failed below even with unlimit stacksize
	pairs = (Pairs*)malloc(sizeof(Pairs)*hold); TEST(pairs)
	n = 0;
	for (i=0; i<secs; i++)
	{ Cells *seci = secs2cell[i]; int si = seci->sort;
		if (si == 2) continue;		// no beta as linker
		for (j=0; j<secs; j++)
	       	{ Cells *secj = secs2cell[j]; int sj = secj->sort; float d;
			if (i <= j) continue;	// no self or double links
			if (secj==seci->sis || secj==seci->bro) continue;
			if (secj->parent != seci->parent) continue; // only intra domain links
			if (si==2 && sj==0) continue;	// no beta-loop
			if (si==0 && sj==2) continue;	// no beta-loop
			if (si==1 && sj==0) continue;	// no alpha-loop
			if (si==0 && sj==1) continue;	// or loop-alpha
			if (si==0 && sj==0) {	// check loop-loop is close
			       if (angle(seci->xyz,seci->parent->xyz,secj->xyz) > PI) continue;
			}
			d = vdif(seci->xyz,secj->xyz);
			if (d > 4.0*bondCA) continue;
			if (d < 2.0*bondCA) continue;
			pairs[n].a = seci;
			pairs[n].b = secj;
			pairs[n].d = d;
			n++;
			if (n==hold) {
				printf("*NB* SSE pairs limit hit\n");
				break;
			}
		}
	}
	qsort(pairs,n,sizeof(Pairs),sortAlpha); // links must be ordered by size to fit distance series
	for (i=0; i<n; i++) { Cells *a = pairs[i].a, *b = pairs[i].b;
		m = a->model*M*N+N; links = data+m+N*5;
		for (j=0; j<links[a->level]; j++) {
			if (a->link[j]==0) { char c,d;
				a->link[j] = b;
				c = d = 'C';
				if (a->sort==1) c = 'A'; if (a->sort==2) c = 'B';
				if (b->sort==1) d = 'A'; if (b->sort==2) d = 'B';
				printf("Linking %d to %d at %f\t(%c%c type in model %d)\n", a->id,b->id,pairs[i].d,c,d,a->model);
				nlinks++;
				break;
			}
		}
	}
	printf("%d links setup\n", nlinks);
}

readlinks ( char *line, FILE *file ) 
{
Cells	*ci, *cj;
FILE	*lin;
int	i,j,k, m,n, in, nlinks, resids = 0, at = 12;
//	read line: LINKS [local|global] <filename> (local = internal numbering, global = PDB resid numbering)
	while (file) {
		in = read_line(file,line);
		if (in < 0) return;
		if (in < 8) continue;
		if (line[0] != 'L') continue;
		if (line[5] != ' ') continue;
		break;
	}
	if (toupper(line[6])=='P') {	// pdbid
		resids = 1;
		printf("Linking by PDB resid numbers\n");
	}
	if (toupper(line[6])=='L') {	// local
		resids = 2;
		printf("Linking by local sequential numbering\n");
	}
	if (toupper(line[6])=='G') {	// global (uid)
		resids = 3; at++;
		printf("Linking by global sequential numbering\n");
	}
	if (resids==0) { Pt(No valid numbering scheme specified) NL return; }
	printf("in: %s\n", line+at);
	lin = fopen(line+at,"r");
	while (1) {
		in = read_line(lin,line);
		if (line[0]=='#') { PRINTs(line+1) NL continue; }
		if (in < 3) { fclose(lin); return; }
		sscanf(line,"%d %d %d", &i, &j, &k);
		if (resids<3) {
			if (resids==1) { i = resn2atom[i]; j = resn2atom[j]; }
			if (resids==2) { i = atom2atom[i]; j = atom2atom[j]; }
			ci = atom2cell[i];
			cj = atom2cell[j];
		} else {
			ci = uid2cell[i];
			cj = uid2cell[j];
		}
		m = ci->model*M*N+N;
		links = data+m+N*5;
		nlinks = links[ci->level];
		if (nlinks>0) { int made = 0;
			for (k=0; k<nlinks; k++) {
				if (ci->link[k] == 0) { // free
					ci->link[k] = cj;
					made = 1;
					break;
				}
			}
			if (made) { float dij = vdif(ci->xyz,cj->xyz);
				printf("Cell %d linked to cell %d   PDBids = %d to %d   dist = %f\n",
					i, j, atom2resn[i], atom2resn[j], dij);
			} else {
				Ps(line) NL
				printf("NB: cell has no free links\n");
			}
		} else {
			Ps(line) NL
			printf("NB: tried to link a cell with no links\n");
		}
	}
}

sceneCell ( Cells *cell )
{
int	i,j,k, n, id;
	id = cell->id;
	n = cell->kids;
	if (n==0) return;
	for (i=0; i<n; i++) {
		for (j=0; j<3; j++) scene->rank[j][total] = cell->child[i];
		total++;
	}
	for (i=0; i<n; i++) sceneCell(cell->child[i]);
	return;
}

modelin ( char *param, int *datat ) {
int	i, j, n, io;
FILE	*dat;
char	*at;
char	line[222];
char	name;
	dat = fopen(param,"r");
	if (!dat) { printf("%s parameter file not found\n", param); exit(1); }
	io = read_line(dat,line);
        if (line[0]=='P') moltype = 0; // protein
        if (line[0]=='R') moltype = 1; // RNA
        if (line[0]=='C') moltype = 2; // Chem (needs fix)
	name = line[io-1];
	datat[0] = moltype;
	// row 0 is internal data not user supplied
	// row 1 is a character tag for each level 
	// specified by an int or all set to name by default
	i = n = 0;
	while (1) {
		if (read_line(dat,line) < 0) break;
		if (line[0] == '/') continue;
		*(strstr(line,"/")) = '\0';
		at = line;
		for (j=0; j<N; j++) { int in;
			sscanf(at,"%d", &in);
			if (in != 0) depth = j;
			datat[N+n] = in;
			printf("%5d ", datat[N+n]);
			n++;
			at = strstr(at,","); at++;
		} NL
		i++;
		if (i==N) break;
	}
	for (i=0; i<N; i++) if(!datat[N+i]) datat[N+i] = (int)name;
	fclose(dat);
	return name;
}

#define T 8

main ( int argc, char** argv ) {
Vec	cog;
Cells	*start;
FILE	*out, *pdb, *run;
char	message[99];
pthread_t thread[T];
int	i, j, m, n, ret[T];
long    rseed = (long)time(0);
        srand48(rseed);
	if (argc>2) { int speed;
		sscanf(argv[2],"%d", &speed);
		if (speed > 0) nanos /= speed;
		if (speed < 0) nanos *= -speed;
		Pt(Running at) Pi(speed) NL
	}
	shrink = 0.0;	// used as 1-shrink (-ve = mass factor)
	scalein = 0.1;	// default unless specified by SCALE <in> <out>
	scaleout = 1.0; // unless specified by SCALE
	data = (int*)malloc(sizeof(int)*4*M*N); TEST(data) // M rows of N numbers for 4 models
	run = fopen(argv[1],"r");
	if (!run) {
		printf("Assume test.run\n");
		run = fopen("test.run","r");
		if (!run) { printf("No test.run file\n"); exit(1); }
	}
	read_line(run,line);
	if (line[0]=='N') {
		delay = -1;		// NORUN flag
		read_line(run,line);
	}
	if (line[0]=='S' && line[1]=='C') { // output SCALE factor (saved in world->far)
		Ps(line) NL
		sscanf(line+7,"%f %f", &scalein, &scaleout);
		read_line(run,line);
	}
	if (line[0]=='S' && line[2]=='R') { // SHRINK factor (saved later in scene->far)
		Ps(line) NL
		sscanf(line+7,"%f", &shrink);
		read_line(run,line);
	}
	if (line[0] != 'P') {
		printf("Need a PARAM filename now\n");
		exit(1);
	}
	nmodels = 0;
	while (1) { char atomtype; // read in param file(s)
		printf("Reading parameters from %s\n", line+6);
		m = nmodels*M*N+N; // read into data array with offset m
		modelin(line+6,data+m);
		// first data row is free (not read in) bar data[m] = moltype (0=PROT, 1=RNA, 2=CHEM)
		names = data+m+N*1; class = data+m+N*2;
        	sizes = data+m+N*3; bumps = data+m+N*4; links = data+m+N*5; chain = data+m+N*6;
		kicks = data+m+N*7; keeps = data+m+N*8; repel = data+m+N*9; bonds = data+m+N*10;
		split = data+m+N*11; // 1 = no link between cousins (set by -bond)
		local = data+m+N*12; // 1 = specific bond length (set by -size)
		align = data+m+N*13; // used by realign()
		atomtype = (char)(names[0]);
		//Pc(atomtype) Pi(depth) NL
		for (i=0; i<N; i++) {
			if (bonds[i]<0) { bonds[i] *= -1; split[i] = 1; } else { split[i] = 0; }
			if (sizes[i]<0) { sizes[i] *= -1; local[i] = 1; } else { local[i] = 0; }
		}
		nmodels++;
		read_line(run,line);
		if (line[0] == 'P') continue;
		if (line[0] == 'M' || line[0] == 'E') { // MODEL or END to finish
			printf("End of parameters\n\n");
			Pt(Read) Pi(nmodels) NL
			for (i=0; i<nmodels; i++) {
				printf("Model %d is", i);
				m = i*M*N+N;
				sizes = data+m+N*3; bonds = data+m+N*10;
				if (data[m]==0) { // protein
					bondCA = 0.1*(sizes[depth]+bonds[depth]);
					Pt(protein with) Pr(bondCA) NL
				}
				if (data[m]==1) { // nucleic
					bondPP = 0.1*(sizes[depth]+bonds[depth]);
					Pt(nucleic with) Pr(bondPP) NL
				}
			}
			break;
		}
		Pt(Bad) Ps(line) NL exit(1);
	}
	if (line[0] == 'M') {
		sscanf(line+6,"%d", &model);
		printf("Using model %d\n", model);
		m = model*M*N+N;
		align = data+m+N*1; class = data+m+N*2;
        	sizes = data+m+N*3; bumps = data+m+N*4; links = data+m+N*5; chain = data+m+N*6;
		kicks = data+m+N*7; keeps = data+m+N*8; repel = data+m+N*9; bonds = data+m+N*10;
		split = data+m+N*11; local = data+m+N*12;
	}
	startup(run);
	fixSSEaxis(world,0,1.0); // fix axes before setangle()
	setangle(world,0,0);	// set the torsion angles and local distances for all levels
	printf("Geom setup done\n");
//exit(1);
	sleep(1);
	/*
	setlinks();      // set intra-chain SSE links
	printf("Link setup done\n");
	atom2pdb = 3.8/(0.1*(float)(sizes[depth]+bonds[depth]));
	atoms = 0;
	out = fopen("given.out","w");
	putpdb(out,world,0,0);
	printf("given PDB saved\n");
	fclose(out);
	// randomise(world,0);
	*/
	presort(data,world,scene);
	printf("Presort done\n");
/*
	for (i=1; i <= depth; i++) {	// fix model for each level to max depth (data[2])
		for (j=0; j<10; j++) {
//			if (align[i]) realign(data,scene,i);
//			n = bondfix(data,scene,i,1);
//			m = bumpfix(data,scene,i,1);
			printf("Level %d: %d bad bonds, %d bad bumps\n", i, n, m);
			if (n+m == 0) break;
		}
		if (j==99) { printf("Input model can't be fixed\n"); exit(1); }
	}
	printf("Patchup done\n");
	setangle(world,0,0);	// reset the torsion angles and local distances for all levels
	if (moltype==0) {
		fixSSEaxis(world,0,1.0);
		fixDOMaxis(world,0,1.0);
		printf("Axis setup done\n");
	}
	vinit(&cog); n = world->kids;
	for (i=0; i<n; i++) vsum(world->child[i]->xyz, &cog);
	vdiv(&cog,(float)n);
	for (i=0; i<n; i++) moveCell(world->child[i],cog,-1); 
	atoms = 0;
	out = fopen("start.out","w");
	putpdb(out,world,0,1);
	printf("start PDB saved\n");
	fclose(out);
	if (delay<0) {
		printf("No run\n");
		chainsout(1); // (1) turns on print
		exit(1);
	} else {
		chainsout(0);
	}
*/
	argcin=argc; argvin=argv;
	focus = (float*)malloc(sizeof(float)*3);
	focus[0] = focus[1] = focus[2] = 0.0;
/* serial debug mode
while (1) {
	acts("Debug acts");
	bump("Debug bump");
	link("Debug link");
	move("Debug move");
}
	sorter(data,world,scene);
	viewer(focus,data,world,scene,types,argcin,argvin);
	sleep(99999);
	sorter(data,world,scene);
	ret[0] = pthread_create( thread+0, NULL, acts, "driver");
	ret[1] = pthread_create( thread+1, NULL, view, "viewer");
	ret[2] = pthread_create( thread+2, NULL, link, "linker");
	ret[3] = pthread_create( thread+3, NULL, move, "mover");
	sleep(99999);
	sorter(data,world,scene);
*/
	// Create independent threads each of which will execute function
/*
for(i=1; i<=22; i++) { Cells *ci = atom2cell[i];
int j; for(j=23; j<=44; j++) { Cells *cj = atom2cell[j];
Pi(i) Pi(j) Pr(vdif(ci->xyz,cj->xyz)) NL
}}
for(i=23; i<=44; i++) { Cells *ci = atom2cell[i];
int j; for(j=1753; j<=1774; j++) { Cells *cj = atom2cell[j];
Pi(i) Pi(j) Pr(vdif(ci->xyz,cj->xyz)) NL
}}
for(i=1; i<=22; i++) { Cells *ci = atom2cell[i];
int j; for(j=3483; j<=3504; j++) { Cells *cj = atom2cell[j];
Pi(i) Pi(j) Pr(vdif(ci->xyz,cj->xyz)) NL
}}
exit(1);
*/
cellout("start.pdb");
	pthread_mutex_init(&linker_lock, NULL);
	pthread_mutex_init(&drawBonds_lock, NULL);
	ret[0] = pthread_create( thread+0, NULL, view, "viewer");
	ret[1] = pthread_create( thread+1, NULL, move, "movers");
/*
	ret[2] = pthread_create( thread+2, NULL, bump, "bumper 0");
	ret[3] = pthread_create( thread+3, NULL, bump, "bumper 1");
	ret[4] = pthread_create( thread+4, NULL, link, "linker");
	ret[5] = pthread_create( thread+5, NULL, rank, "sorter");
	ret[6] = pthread_create( thread+6, NULL, fixs, "fixers");
	ret[7] = pthread_create( thread+7, NULL, acts, "driver"); // NB acts() counts the frames 
*/
	sleep(99999);
	pthread_mutex_destroy(&drawBonds_lock);
	pthread_mutex_destroy(&linker_lock);
/*
	m = 0;
	while (1) { Vec cog; 
		sleep(1);
		setangle(world,0,1);		// update local ideal angles and distances
		sleep(1);
		for (i=1; i <= depth; i++) {	// for each level to max depth (data[2])
			if (align[i]) realign(data,scene,i);
			sleep(1);
			bondfix(data,scene,i,0);
			sleep(1);
			bumper(data,world,scene,i);
		}
		sleep(1);
		fixSSEaxis(world,0,0.5);
		sleep(1);
		fixDOMaxis(world,0,0.5);
		sleep(1);
		vinit(&cog); n = world->kids;
		for (i=0; i<n; i++) vsum(world->child[i]->xyz, &cog);
		vdiv(&cog,(float)n);
		for (i=0; i<n; i++) moveCell(world->child[i],cog,-1); 
		sleep(1);
		chainsout(1);
		sleep(1);
		atoms = 0;
		out = fopen("test.out","w");
		putpdb(out,world,0,1);
		printf("PDB saved\n");
		sleep(delay);
		m++;
		if (m==10) break;
		if (m==50) break;
		kicks = data+N+N*7;
		kicks[2]/=2; kicks[3]/=2;
		Pi(kicks[1]) Pi(kicks[2]) Pi(kicks[3]) NL
		repel = data+N+N*9;
		repel[1]/=2; repel[2]/=2;
		Pi(repel[1]) Pi(repel[2]) Pi(repel[3]) NL
	}
*/
	printf("Exiting\n");
	//for (i=0; i<5; i++) pthread_join(thread[i], NULL);
}

void cellout ( char *filename ) {
int	i, n, lev, dom;
float	fsort;
FILE	*out;
char	c;
Pi(total) Pr(scaleout) NL
	out = fopen(filename,"w");
	DO(i,total) { Cells *cell = uid2cell[i];
		if (cell->level==0) continue;
		lev = depth - cell->level;
		if (lev) fsort = cell->sort; else fsort = cell->parent->sort;
		if (fsort > 3.0) fsort = -1.0;
		if (cell->level < 2) dom = -1; // act/myo
		if (cell->level==2) dom++;
		if (dom<0) c = 'A'; else c = 'A'+dom;
		if (dom>26) c = 'a'+ dom-26;
		n = 0;
		if (cell->level==depth) n = cell->resn;
		if (n>9999) n = 9999;
fsort = log((float)(cell->busy+1));
		if (lev) fprintf(out,"ATOM%7d  C%cn GLY %c%4d     %7.3f %7.3f %7.3f %5.2f %5.2f\n", i, 'A'+lev, c, n,
			scaleout*cell->endN.x, scaleout*cell->endN.y, scaleout*cell->endN.z, (float)cell->type, fsort);
        	fprintf(out,"ATOM%7d  C%c  GLY %c%4d     %7.3f %7.3f %7.3f %5.2f %5.2f\n", i, 'A'+lev, c, n,
			scaleout*cell->xyz.x, scaleout*cell->xyz.y, scaleout*cell->xyz.z, (float)cell->type, fsort);
        	if (lev) fprintf(out,"ATOM%7d  C%cc GLY %c%4d     %7.3f %7.3f %7.3f %5.2f %5.2f\n", i, 'A'+lev, c, n,
			scaleout*cell->endC.x, scaleout*cell->endC.y, scaleout*cell->endC.z, (float)cell->type, fsort);
		if (cell->level==2) fprintf(out,"TER\n");
	}
}

void chainout ( int print, FILE *out, FILE *pdb, Cells *cell, int id)
{
int	i, m, rna, forks, end=0, level = cell->level;
float	type, sort, sortout, scale = 0.3;
	if (cell->done) return;
	m = cell->model*M*N+N;
	if (data[m]) rna = 1; else rna = 0;
        sizes = data+m+N*3; bumps = data+m+N*4; links = data+m+N*5; chain = data+m+N*6;
	kicks = data+m+N*7; keeps = data+m+N*8; repel = data+m+N*9; bonds = data+m+N*10;
	forks = abs(chain[cell->level]);
	type = (float)cell->type;
	sort = (float)cell->sort;
	if (cell->parent->type == 2) {
		type = (float)cell->parent->type;
		sort = (float)cell->parent->sort;
	}
	atoms++;
	if (sort < 2.5) sortout = 1.0;
	if (sort < 1.5) sortout = 2.0;
	if (sort < 0.5) sortout = 0.0;
        fprintf(out,"ATOM%7d  CA  GLY %c%4d     %7.3f %7.3f %7.3f %5.2f %5.2f\n", atoms, 'A'+id, atoms,
		scale*cell->xyz.x, scale*cell->xyz.y, scale*cell->xyz.z, type, sortout);
	if (cell->level == depth) { float x,y,z, s = atom2pdb, t = 400;
		if (rna) s *= 0.65; // scale RNA to fit fake protein in SAP
		x = t+s*cell->xyz.x; y = t+s*cell->xyz.y; z = t+s*cell->xyz.z;
		if (x < -99.999) x = -99.999; if (x > 999.999) x = 999.999;
		if (y < -99.999) y = -99.999; if (y > 999.999) y = 999.999;
		if (z < -99.999) z = -99.999; if (z > 999.999) z = 999.999;
		fprintf(pdb,"ATOM%7d  CA  GLY %c%4d     %7.3f %7.3f %7.3f %5.2f %5.2f\n",
			atoms, 'A'+i, atoms, x, y, z, type, sortout);
	}
	if (print && cell->bond)
	{ Cells	*csis = cell->sis, *cbro = cell->bro;
	  float	dn = vdif(cell->xyz,csis->xyz);
	  float	dc = vdif(cell->xyz,cbro->xyz);
	  float danger = 2.0*0.1*(float)(sizes[level]+bonds[level]);
		if (cell == cell->parent->starts) end += 1;
		if (cell == cell->parent->finish) end += 2;
		printf("%4d  %1d  %4d%4d%4d  %4d%4d%4d  %4d%4d%4d ", level,end,
			cell->atom,cell->id,cell->parent->id,
			csis->atom,csis->id,csis->level,
			cbro->atom,cbro->id,cbro->level);
		if (level==csis->level) Pr(dn)
		if (level==cbro->level) Pr(dc)
		if (dn>danger || dc>danger) printf(" ********");
		NL
	}
/*
	cell->done = 1;	// set output mark
	for (i=0; i<forks && cell->bond && cell->bond[i][1]; i++) 
	{ Cells *next = cell->bond[i][1];
		if (next->level != level) continue;
		chainout(print, out, pdb, next, id);
	}
*/
}
/*
findStart ( int print, FILE *out, FILE *pdb, Cells *cell, int level )
{
int	i,j,k, m,n, id, forks;
	n = cell->kids;
	if (level && cell->level==level && !cell->done) { // at right level and not done
		m = cell->model*M*N+N; chain = data+m+N*6;
		forks = abs(chain[level]);
		if (print) printf("Level %d chain %d\n", level, cell->id);
		if (cell->bond) { 	// in a chain
			for (i=0; i<forks; i++) 
			{ Cells *bond = cell->bond[i][0];
				if (bond && bond->level<level) {
					if (print) printf("level %d chain %d\n", level, cell->id);
					atoms = 0;
					chainout(print,out,pdb,cell,level-1);
					if (level==depth) fprintf(pdb,"TER\n");
					fprintf(out,"TER\n");
					break;
				}
			}
		} else {		// isolated cell
		}
	}
	for (i=0; i<n; i++) findStart(print,out,pdb,cell->child[i],level);
	return;
}

chainsout (int print)
{
FILE	*out, *pdb;
int	level;
	out = fopen("chain.out","w");
	pdb = fopen("chain.pdb","w");
	for (level=1; level<=depth; level++) {
		findStart(print,out,pdb,world,level); // check for free starts
	}
	printf("PDB chains saved\n");
	fclose(out);
	fclose(pdb);
}
*/
 
void putpdb ( FILE* out, Cells* cell, int level, int ends )
{
float	type, sort, scale = 3.8/(0.1*(float)(sizes[depth]+bonds[depth]));
float	sortout, d, err, max, bond, bondlen, bonderr;
int     i, bad, n = cell->kids;
int	m = cell->model*M*N+N;
int	rna = data[m];
       	sizes = data+m+N*3; bumps = data+m+N*4; links = data+m+N*5; chain = data+m+N*6;
	kicks = data+m+N*7; keeps = data+m+N*8; repel = data+m+N*9; bonds = data+m+N*10;
	if (moltype==1) scale *= 0.2; // smaller for RNA (so points join in RASMOL)
	cell->done = 0;
	if (n==0) return;
	if (level == depth-1) { // write child atoms only at lowest level
		for (i=0; i<n; i++) { Cells *child = cell->child[i];
			type = (float)cell->type;
			sort = (float)cell->sort;
			if (sort < 2.5) sortout = 1.0;
			if (sort < 1.5) sortout = 2.0;
			if (sort < 0.5) sortout = 0.0;
			atoms++;
                	fprintf(out,"ATOM%7d  CA  GLY%6d     %7.3f %7.3f %7.3f %5.2f %5.2f\n", atoms, atoms,
				scale*child->xyz.x, scale*child->xyz.y, scale*child->xyz.z, type, sortout);
		}
		// write parent position
               	fprintf(out,"ATOM%7d  CB  GLY%6d     %7.3f %7.3f %7.3f %5.2f %5.2f\n", atoms, atoms,
			scale*cell->xyz.x, scale*cell->xyz.y, scale*cell->xyz.z, type, sortout);
		if (ends && sort>0.5 && type>1.5) { // write ends if in sec and not sphere
               		fprintf(out,"ATOM%7d  N   GLY%6d     %7.3f %7.3f %7.3f %5.2f %5.2f\n", atoms, atoms,
				scale*cell->endN.x, scale*cell->endN.y, scale*cell->endN.z, type, sortout);
               		fprintf(out,"ATOM%7d  O   GLY%6d     %7.3f %7.3f %7.3f %5.2f %5.2f\n", atoms, atoms,
				scale*cell->endC.x, scale*cell->endC.y, scale*cell->endC.z, type, sortout);
		}
	}
	if (cell->bond && chain[level+1]) {
		bondlen = bonderr = max = 0.0;
		bond = 0.1*(float)(sizes[level+1]+bonds[level+1]);
		for (i=1; i<n; i++) { Cells *child = cell->child[i];
			d = vdif(child->xyz,child->sis->xyz);
			bondlen += d;
			d = d - bond; 
			err = d*d;
			if (err > max) { max = err; bad = i; }
			bonderr += err;
		}
		bondlen /= (float)(n-1);
		bonderr /= (float)(n-1);
		bonderr = sqrt(bonderr);
		if (level < depth-1 && bonderr > NOISE*10.0) {
			d = vdif(cell->child[bad]->xyz, cell->child[bad]->sis->xyz);
///			printf("Worst bond length for level %d bond %d--%d = %f for ideal %f (mean %f)\n",
///				level+1, bad-1, bad, d, bond, bondlen);
		}
	}
	for (i=0; i<n; i++) putpdb(out, cell->child[i], level+1, ends);
}
 
void fixRNAstem ( Cells* cell, int level, float weight )
{
Vec	*new1, *new2;
Vec	p,q,r, u,v,w, at, dt, axis, shift;
int     i, j, m, n = cell->kids;
float	rise, radius, twist, h, ab,bc,cd, ad;
float	clash = 1.6;
	m = cell->model*M*N+N; sizes = data+m+N*3;
	moltype = data[m];
	radius = 0.05*(float)sizes[level];
	if (n==0) return;
	if (moltype==1 && cell->type==2 && cell->sort) { // RNA stem
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
/*
if (i) {
float dis = vdif(a,cell->child[i-1]->xyz),
bond = 0.1*(float)(sizes[level+1]+bonds[level+1]);
bc = vdif(b,c);
ad = vdif(a,d);
Pr(dis) Pr(bond) Pr(radius) Pr(ab) Pr(cd) Pr(bc) Pr(ad) NL
}
*/
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
			DO(i,n) {
				cell->child[i]->lost = 1;
				mids[i].x = 999.9;
			}
			DO(i,n) { Cells *ci = cell->child[i];
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
			DO(i,k) { float d, dmin;
				// shift newN towards nearest remaining mid-point
				dmin = 999.9;
				DO(j,n) {
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
				DO(j,n) {
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
			DO(j,bad+bent) repair(cell);
		}
		DO(i,n) cell->child[i]->lost = 0;
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
/*
Pi(mid) NL
i = 0;
printf("ATOM%7d  CA  ARG A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,cell->endN.x,cell->endN.y,cell->endN.z);
i = 1;
printf("ATOM%7d  CA  GLU A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,cell->endC.x,cell->endC.y,cell->endC.z);
printf("TER ATOM\n");
DO(i,len) { Cells *ci = cell->child[i];
printf("ATOM%7d  CA  CYS A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,ci->id,ci->xyz.x,ci->xyz.y,ci->xyz.z);
}
printf("TER ATOM\n");
*/
//rad = 0.76;
	vsub(cell->endC,cell->endN,&x); vnorm(&x);
	vcopy(x,&r); vmul(&r,turn/11.0); // rise/base
	vcopy(x,&t); vmul(&t,turn*0.5); // half turn
	DO(j,len) {
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
	DO(jj,n) { int new;
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
DO(i,len) { Cells *ci = cell->child[i];
printf("ATOM%7d  CA  CYS A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,ci->xyz.x,ci->xyz.y,ci->xyz.z);
}
printf("TER ATOM\n");
	vsub(cell->endC,cell->endN,&x); vnorm(&x);
	vcopy(x,&r); vmul(&r,rise);
	vcopy(x,&t); vmul(&t,turn);
	DO(j,len)
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
/*
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
*/
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
DO(i,len) { Cells *ci = cell->child[i];
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
	DO(i,n) { Cells *ci = cell->child[i];
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
 
void fixSSEaxis ( Cells* cell, int level, float weight )
{
Vec	x, y, this, last, axis, shift;
int     i, k, m, in, n = cell->kids;
float	fn = (float)n;
	m = cell->model*M*N+N; sizes = data+m+N*3;
	moltype = data[m];
	if (n==0) return;
	if (moltype==0 && cell->type==2) { // protein SSE (and not RNA)
		if (cell->sort) { // alpha or beta (set line by i-1,i+1 pairs)
			in = 0;
			vinit(&axis);
			vave(cell->child[2]->xyz,cell->child[0]->xyz, &last);
			for (i=2; i<n-1; i++) {
				vave(cell->child[i-1]->xyz,cell->child[i+1]->xyz, &this);
				vsub(this,last, &shift);
				vsum(shift, &axis);
				vcopy(this,&last);
				in++;
			}
			if (in<1) return;
			vnorm(&axis);
			if (cell->sort==1) vmul(&axis,0.5*alphAXIS*fn*bondCA/3.8);
			if (cell->sort==2) vmul(&axis,0.5*betaAXIS*fn*bondCA/3.8);
		} else { // loop
			if (n<2) return;
			m = n/2;
			vinit(&x); k = 0;
			for (i=0; i<m; i++) {
				vsum(cell->child[i]->xyz,&x); k++;
			}
			vdiv(&x,(float)k);
			vinit(&y); k = 0;
			for (i=m; i<n; i++) {
				vsum(cell->child[i]->xyz,&y); k++;
			}
			vdiv(&y,(float)k);
			vsub(y,x, &axis);
			vnorm(&axis); vmul(&axis,0.5*loopAXIS*sqrt(fn));
		}
		vsub(cell->xyz, axis, &x);
		vsub(x,cell->endN, &shift); vmul(&shift,weight);
		vsum(shift, &(cell->endN));
		vadd(cell->xyz, axis, &y);
		vsub(y,cell->endC, &shift); vmul(&shift,weight);
		vsum(shift, &(cell->endC));
		cell->ends = 1;
/*
for (i=0; i<cell->kids; i++) { Cells *ci = cell->child[i];
printf("ATOM%7d  CA  GLY A%4d     %7.3f %7.3f %7.3f  0.00 %5.2f\n", i,i,ci->xyz.x,ci->xyz.y,ci->xyz.z,(float)ci->lost);
} i=0;
printf("TER ATOM\n");
printf("ATOM%7d  CA  CYS A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,cell->endN.x,cell->endN.y,cell->endN.z);
printf("ATOM%7d  CA  CYS A%4d     %7.3f %7.3f %7.3f  0.00  0.00\n", i,i,cell->endC.x,cell->endC.y,cell->endC.z);
printf("TER ATOM\n");
if(cell->id==2) exit(1);
*/
	}
	for (i=0; i<n; i++) fixSSEaxis(cell->child[i], level+1, weight);
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
/*
	} else {
		if (vdif(cell->endN,cell->endC) > NOISE) { // endpoints exist
			cell->ends = 1;
		}	
*/
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
		DO(i,n)
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
	DO(i,n) reset2axis(cell->child[i], level+1, weight);
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
				DO(i,n) { Cells *kidi = cell->child[i]; // set cent reference distances
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
		DO(i,n)
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
	DO(i,n) fixANYaxis(cell->child[i], level+1, weight);
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
