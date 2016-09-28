#define M 5	// maximum number of models
#define N 10	// size of model data array (N cols per row = levels))
#define E 20    // sorts of ellipsoids
#define SPARE 3 // extra work/scratch copies off the end of the world
// protein SSE dimensions: step along AXIS/residue, relative THICkness
#define alphAXIS 1.5
#define betaAXIS 3.0
#define loopAXIS 0.5
#define alphTHIC 1.2
#define betaTHIC 0.8
#define loopTHIC 1.0
#define CALM 2 // used to dec. busy in sorter() and cell->bump in viewer()

class Data {
public:
	int	moltype;	//  molecule type
	int	subtype;	//  molecule subtype 
	int	colours;	//  colour scheme
	// model marameters
	char	names[N];	//  not used
	int	shape[N];	//  shape type: 0=virtual, 1=sphere, 2=tube, 3=ellipsoid
				//	-ve = render as wireframe
	float	sizes[N];	//  size that appears when rendered and used by keeper()
				//	-ve = use values in prox and dist
	float	bumps[N];	//  size used by bumper()
				//	-ve
	int	links[N];	//  number of link[] slots allocated
				//	-ve
	int	chain[N];	//  topology: number of bonds (+ve = linear)
				//	-ve = circular
	float	kicks[N];	//  kick strength of a move (used as kicks[]*0.001)
				//	-ve = no spin
	float	keeps[N];	//  keep strength (used as keeps[]*0.01)
				//	-ve = shell for types 1,3 or tube = open 
	float	bonds[N];	//  bond length between objects' surface (size)
				//	-ve
	float	repel[N];	//  kick strength of a bump repulsion
				//	-ve don't repel children within parent 
	float	rejel[N];	//  kick strength for soft repulsion (jelly)
				//	-ve don't repel children within parent 
	int	split[N];	//  flag to mark no bond between cousins (set by -ve bond length)
	int	local[N];	//  flag for generic=0 or specific=1 bond length (set by -ve size)
	// global values
	static	Data	model[M];	// the set of Data models
	static	float	focus[3];	// focus point of the view
	static	float	Eratio[E];	// range of ellipsoid shapes
	static	float	shrink;
	static	float	bondCA, bondPP;
	static	float	scalein, scaleout;
	static	int	depth, nmodels;
	static	int	norun, noview;
	static	int	frame, hidden;
};

extern char  *names;
extern int   *shape, *links, *chain, *split, *local, *align;
extern float *sizes, *bumps, *kicks, *keeps, *bonds, *repel, *rejel;
extern int   total, depth, model, moltype, subtype;
