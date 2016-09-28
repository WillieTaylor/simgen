#include "util.hpp"
#include "geom.hpp"
#include "cell.hpp"
#include "data.hpp"

int argcin; char** argvin;

void params ( FILE* );
void models ( FILE* );
void bumper ();
void shaker ();
void potter ( int, int );
void tinker ( float, int, int, int );
void keeper ( Cell* );
void bonder ( Cell* );
void linker ( Cell* );
void center ();
void sorter ();
void viewer ( int, char** );
void presort ();
void driver (); // main application code
void helper (); // application code addition
void minder (); // application code patchups
void looker (); // application code monitors
void setter (); // application code initialise
Cell* setScene ( Cell* );
void fixSSEaxis( Cell*, float );
void fixRNAaxis( Cell*, float );
void fixRNAstem( Cell*, float );

//           100000000 = 1/10 sec.
long nanos = 100000000;	// good speed for watching
int  delay = 1;//10;
int  lastlook = 0;

float	atom2pdb;	// scale factor

void *acts ( void *ptr )
{
	sleep(1);	// delay start to allow stuff to be set
	setter();	// one-off application init (user provided)
	sleep(1);	// delay start to allow stuff to be set
	printf("Starting driver\n");
	while (1) {
		driver(); // application control (user provided)
		Cell::world->idata[0]++;
	}
}

void *look ( void *ptr )
{
	sleep(1);
	printf("Starting looker \n");
	while (1) {	// looker() should contain its own timeout delay
		looker(); // watches and measures things (user provided)
		helper(); // helps the driver with stuff (user provided)
		lastlook++;
	}
}

void *fixs ( void *ptr )
{
float	weight = 0.1; //(float)nanos/100000000.0; // = 1 when speed = 1, .1 for 10, .01 for 100
	sleep(1); // delay start to allow ends to be set
	printf("Starting fixers \n");
	while (1) {	// minder() should contain its own timeout delay
		minder(); // helps the driver with stuff (user provided)
//		fixSSEaxis(Cell::world,weight);
		weight = Data::fixRNA;	// can be changed in used code (eg: minder())
		fixRNAstem(Cell::world,weight);
		timeout();
		fixRNAaxis(Cell::world,weight);
/*
		nanosleep(&wait,NULL);
		fixRNAstem(world,0,weight);
		nanosleep(&wait,NULL);
		fixDOMaxis(world,0,1.0);
		fixANYaxis(world,0,weight);
		nanosleep(&wait,NULL);
		reset2axis(world,0,weight);
*/
		Cell::world->idata[1]++;
		if (Data::frozen) // update the focus position and heat values
		{ Cell	*focus = Cell::uid2cell[abs(Data::frozen)];
		  float range = 0.1, size = 0.5*focus->len; size *= size;
			Data::focus = focus->xyz;
			FOR(i,Cell::total+1) // set heat from the squared distance from the focus
			{ Cell	*c = Cell::uid2cell[i]; float dd;
				if (c->heat > 1.0) continue; // flag for exempt level (over focus)
	  	  		dd  = c->xyz || Data::focus;
				if (dd < size) {
					c->heat = 1.0;
				} else {
					c->heat = exp(-dd*range/size);
				}
			}
		}
	}
}

void *move ( void *ptr )
{
	printf("Starting mover\n");
	while (1) { int fix = abs(Data::tinker);
		Data::frame++;
		shaker();		// makes random displacements at all levels
		keeper(Cell::world);	// keeps children inside their parent
		center();		// keeps children centred on their parent
		if (fix>0 && (Data::frame%fix)==0) {
			tinker(0.5,Data::depth,0,1);	// geometry patch-up using stored local data
		}
		Cell::world->idata[2]++;
		timeout(nanos);
	}
}

void *link ( void *ptr )
{
	printf("Starting linker \n");
	while (1) {
		bonder(Cell::world);
		linker(Cell::world);
		Cell::world->idata[3]++;
		timeout(nanos);
	}
}

void *rank ( void *ptr )
{
	printf("Starting sorter \n");
	while (1) {
		sorter(); // update the ranked XYZ lists
		//timeout(nanos);
		timeout(1,nanos); // <-------NB + 1 sec.
	}
}

void *bump ( void *ptr )
{
	printf("Starting bumper \n");
	while (1) {
		bumper();
		Cell::world->idata[4]++;
		timeout(nanos/10);
	}
}

void *pots ( void *ptr )
{
struct  timespec start, finish;
long	runtime, delay;
int	m = Data::depth;
char	*flag = (char*)ptr;
int	cycles = 1;
	printf("Starting molecular dynamics for %s\n", flag);
	while(1) { timespec start,finish;
   		clock_gettime(CLOCK_REALTIME, &start);
		if (flag[0]=='a') {
			potter(m,cycles);
		} else {
			DO1(i,m-1) potter(i,(cycles*m)/i);
		}
		clock_gettime(CLOCK_REALTIME, &finish);
		runtime = timedifL(start,finish);
		if (runtime < nanos) {
			delay = nanos - runtime;
			//Ps(flag) Pt(waiting for) Pi(delay) Pi(cycles) NL
			timeout(delay);
			cycles++;
		}
	}
}

void *view ( void *ptr )
{
	printf("Starting viewer \n");
	viewer(argcin,argvin); // sets-up graphical objects (types) then runs by callback
}

main ( int argc, char** argv ) {
pthread_t thread[10];
FILE    *run;
int	i, j, m, n;
long    rseed = (long)time(0);
        srand48(rseed);
	if (argc>2) { int speed;
		sscanf(argv[2],"%d", &speed);
		if (speed > 0) { nanos /=  speed; Data::delay /=  speed; }
		if (speed < 0) { nanos *= -speed; Data::delay *= -speed; }
		Pt(Running at) Pi(speed) Pr(Data::delay) NL
	}
	run = 0;
        if (argc > 1) run = fopen(argv[1],"r");
        if (!run) {
                printf("Assume test.run\n");
                run = fopen("test.run","r");
                if (!run) { printf("No test.run file\n"); return 1; }
        }
	Cell::world = new Cell;		// create the top-level cell structure (set in models())
	params(run);
	models(run);
	fclose(run);
	Pi(Cell::total) NL
	Pi(Data::depth) NL
	Pr(Data::shrink) NL
	Cell::world->print();
	Cell::world->move(-Cell::world->xyz);
	Pt(Reset) Pv(Cell::world->xyz) NL
	Cell::scene = setScene(Cell::world);
	Cell::scene->sort = 0;
	presort();
	Data::frame = -10;
	argcin=argc; argvin=argv;
	if (Data::noview) {
		Pt(No viewing) NL
	} else {
		pthread_create( thread+0, 0, view, 0 );
	}
	pthread_create( thread+1, 0, rank, 0 );
	if (Data::nomove==0) { // skipped on NOMOVE
		pthread_create( thread+2, 0, move, 0 ); // shaker
		pthread_create( thread+3, 0, bump, 0 ); // bumper
		pthread_create( thread+4, 0, link, 0 ); // linker and bonder
		pthread_create( thread+5, 0, fixs, 0 ); // fixers and minder
	}
	if (Data::nolive==0) { // skipped when NORUN set
		pthread_create( thread+6, 0, acts, 0 ); // driver
		pthread_create( thread+7, 0, look, 0 ); // looker and helper
	}
/* toy MD calls
	pthread_create( thread+8, 0, pots, (void*)"atoms" );
	pthread_create( thread+9, 0, pots, (void*)"other" );
*/
	sleep(99999);
}

/* code to stretch chain
if (Data::frame<200) {
FOR(i,Cell::total) { Cell *c=Cell::uid2cell[i];
if (c->level==3) c->xyz.x -= 0.005*(i-Cell::total/2);
}
FOR(i,Cell::total) { Cell *c=Cell::uid2cell[i];
if (c->level==2) { Vec x = (c->endC - c->endN)*0.5;
c->xyz = c->starts->xyz & c->finish->xyz;
c->endN = c->xyz-x; c->endC = c->xyz+x; }
} 
nanosleep(&wait,NULL);
continue;
}
*/
