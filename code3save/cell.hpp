#define NDATA 5 
#define HOLD 10        // HOLD x kids = maximum number of potential bumping pairs held
#define LIVE 10 // 1

struct bonds;

class Cell {
public:
        Vec     xyz;		// coordinates
	int	model;		// model type
	Cell	*parent;	// upper level
	Cell	*sis, *bro;	// adjacent siblings
	Cell	**child;	// list of children
	Cell	*junior;	// youngest child [0]
	Cell	*senior;	// eldest child [n-1]
	Cell	*starts;	// child containing chain start (can = finish)
	Cell	*finish;	// child containing chain finish (can = start)
	Cell	***rank;	// list of members ranked in XYZ (for each dimension)
	struct	bonds	*link;		// list of links
	struct	bonds	*bond;		// list of bonds (pointers to cells)
	int	nbonds, nlinks, cots;	// space (re)allocated for kids
	int	ranks[3];		// position in parents XYZ list
	int     id, kids, uid, colour;	// id = number in family kids, uid = unique id
	int	level, type, sort, ends;// sort of the type found at level (with ends)
	Vec	endN, endC, geom;	// connection points, local angles
	Vec	prox, dist, cent;	// local distances (proximal, distal, centroid)
	Mat	rot;			// rotation frame for body
	int	atom, btom, ctom, resn;	// low-level, sub-level, all atom number
	int	idata[NDATA];		// int flags/values for driver to use
	float	fdata[NDATA];		// float flags/values for driver to use
	Vec	vdata[NDATA];		// vector flags/values for driver to use
	char	cdata[NDATA];		// char flags/values for driver to use
	float	far, len, heat;		// sq-dist to eye, axis length, sq-dist to focus
	int     solid, lost, empty;	// flags for rendering
	int	nsis, nbro;		// bond counters
	int	live, done, bump, busy; // flags for active/bumping cell
	Cell	*hit;		// bumping cell
	// global data
	static	Cell	*world;
	static	Cell	*scene;
	static	int	total, allatoms;
	static	Cell	**uid2cell;
	// constructors
	Cell() {
		id = kids = level = empty = solid = ranks[0]=ranks[1]=ranks[2] = 0;
		live = 1;
	}
	// functions
	void extend ( const int, const int, const int );	// extend the size of the current cell
	void spawn ( const int, const int, const int );	// create n child objects for the current cell
	void spawn ( const int, const int );	// create n child objects for the current cell
	void spawn ( const int );		// create n child spheres for the current cell
	void spawn ();				// create one sphere for the current cell
	void print ();				// print the cell and decendrants
	void move ();				// move the cell and decendants in a random unit direction
	void move ( const float );		// move the cell and decendants in a random direction by d
	void move ( const Vec );		// move the cell and decendants
	void rankXYZ ( const int );		// rank the XYZ lists of child positions
	int  bumpin ();				// bumping children within a cell
	void bumps  ();				// bumps within current and all offspring
	void rankXYZ () {  rankXYZ(10); }	// default 10 sort passes of rankXYZ()
	void sortXYZcell ( const int );		// pass over cell tree calling sortXYZrank()
	int  sortXYZrank ( const int );		// one pass bubble sort on the children of the current cell
	Vec  getWcent ();			// returns the weighted centroid
	void setWcent ();			// sets the weighted centroid
	void group();				// lock children centroid to parent centre
	void shifter( const float );		// couples parent's centre to centroid of kids
	void spin ( const Vec zero, const Vec axis, const float theta ); // spin cell tree 
	void spin ( const Vec axis, const float theta ); // spin about centre--axis
	void spin ( const Seg line, const float theta ); // spin about a line
	void spin ( const float theta );		// spin about a random line
	void spin ();				// spin about a random line and angle
};

struct bonds {
	struct	Cell	*to;	// bonded cell
	int	type;		// bond thickness(0--3 for drawLink())
	int	next;		// next bond to follow from "to"
	int	link;		// connect 0=1=NN, 2=NC, 3=CN, 4=CC
}; typedef struct bonds Bonds;	// indirect typedef because bonds is used in Cell

typedef struct {
        Cell   *a, *b;
        int     c;
        float   d;
        float   bump;
} Bumps;

///////// utility

void moveCell ( Cell*, Vec, int );	// same as Cell::move() but keeps old style parameter list
void moveCell ( Cell*, Vec );		// use positive shift

void spinCell ( Cell*, Vec, Vec, float );	// same as Cell::spin() but old style parameter list
void turnCell ( Cell*, Vec, Mat );	// rotates cell about point vec using mat

void part2cells ( Cell*, Cell*, float, float ); // separates two cells towards target distance
void part2cells ( Cell*, Cell*, float );	// use full separation (kick = 1)

void setWcent ( Cell* );		// sets the weighted centroid in cell
Vec  getWcent ( Cell* );		// returns the weighted centroid of cell

///////// bumper

bool exempt ( Cell*, Cell* );		// TRUE when childern are exempt from bumping
int bumpex ( Cell*, Cell* );		// bumps children between two cells
float touch ( Cell*, Cell* );		// closest approach between two cell surfaces
float touch ( Cell*, Cell*, int, int );	// closest approach between two cell (bump/shell) surfaces
int getBumpin ( Cell*, Bumps* );	// return the possible bumping pairs in the cell

///////// check for links

int  getLink  ( Cell*, Cell* );		// returns the link(+1) from a to b (-ve b to a)
bool isLinked ( Cell*, Cell* );		// check if any link between a and b
bool kidsLink ( Cell*, Cell* );		// check if any kids link between A and B
bool kidsLink ( Cell*, int, Cell* );	// check if kid a in A has a link with B
bool kidsLink ( Cell*, int, Cell*, int );	// check if kid Aa links with Bb

float tube_to_egg ( Cell*, Cell* );
float  egg_to_egg ( Cell*, Cell* );

///////// print (in pdb format)
				// all levels unsorted with scale
void dumpall ( float );			// with scale
				// all levels (sorted as separate chains)
void sortpdb ( float );			// with scale and pdbid numbering
void sortall ( float );			// with scale
void sortall ();				// default: scaleout
				// just atomic level
void putpdb ( FILE*, Cell*, float );	// into open file
void putpdb ( char*, Cell*, float );	// open new filename
void putpdb ( Cell* );			// default "dump.pdb", scaleout
void putpdb ( Cell*, float );		// default "dump.pdb"
void putpdb ( char*, float );		// default: world
void putpdb ( char* );			// default: world, scaleout
				// print chains in (re)linked order
void pdbout ( FILE*, Cell*, float );	// into open file
void pdbout ( char*, Cell*, float );	// open new filename
void pdbout ( char*, Cell* );		// dito using default scaleout
void pdbout ( char* );			// dito using world
