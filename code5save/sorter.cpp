#include "util.hpp"
#include "geom.hpp"
#include "cell.hpp"
#include "data.hpp"

typedef struct {
        Cell   *a, *b;
        int     c;
        float   d;
} cellPairs;

int sortPairs ( const void *ac, const void *bc )
{
        cellPairs   *a = (cellPairs*)ac, *b = (cellPairs*)bc;
        if (a->d < b->d) return -1;
        if (a->d > b->d) return  1;
        return 0;
}

int sortXYZview ( int, int );

void sorter ()
{
int	i, j, n, resort = 10;
int	use = 2, end = 0;
int	sortcode = Cell::scene->sort;
	// sort internal lists
	for (j=0; j<10; j++) {
		for (i=0; i<3; i++) Cell::world->sortXYZcell(i);
	}
	// sort all view lists
	for (i=0; i<3; i++) {
		for (j=0; j<resort; j++) { if (!sortXYZview(i, 1)) break; }
		for (j=0; j<resort; j++) { if (!sortXYZview(i,-1)) break; }
	}
	// sort in view direction (data[5] : 123=XYZ, +/-=head/tail)
	for (j=0; j<resort; j++) { int use = 2, view = 1;
		if (sortcode > 0) { use =  sortcode-1; view =  1; }
		if (sortcode < 0) { use = -sortcode-1; view = -1; }
		if (!sortXYZview(use,view)) break;
	}
}

int sortXYZview (int use, int flip)
{
Cell	*scene = Cell::scene;
int	total = Cell::total;
int	i, swap = 0;
float	ga,gb, ra,rb;
	if (flip > 0) {	// shift forward by a radius for front view
		for (i=1; i<total; i++)
		{ Cell	*a = scene->rank[use][i-1], *b = scene->rank[use][i];
		  float sizea = Data::model[a->model].sizes[a->level],
		  	sizeb = Data::model[b->model].sizes[b->level];
			// returning 0 stops sorting cycle
			if (a==0 || b==0) return 0;
			if (use==0) { ga = a->xyz.x; gb = b->xyz.x; }
			if (use==1) { ga = a->xyz.y; gb = b->xyz.y; }
			if (use==2) { ga = a->xyz.z; gb = b->xyz.z; }
			ga += sizea; gb += sizeb;
			if (ga > gb) {
				scene->rank[use][i] = a;
				scene->rank[use][i-1]=b;
				swap++;
			}
		}
	} else {	// shift backward by a radius for rear view
		for (i=total+1; i<total*2; i++)
		{ Cell *a = scene->rank[use][i-1], *b = scene->rank[use][i];
		  float sizea = Data::model[a->model].sizes[a->level],
		  	sizeb = Data::model[b->model].sizes[b->level];
			// returning 0 stops sorting cycle
			if (a==0 || b==0) return 0;
			if (use==0) { ga = a->xyz.x; gb = b->xyz.x; }
			if (use==1) { ga = a->xyz.y; gb = b->xyz.y; }
			if (use==2) { ga = a->xyz.z; gb = b->xyz.z; }
			ga -= sizea; gb -= sizeb;
			if (ga < gb) {
				scene->rank[use][i] = a; 
				scene->rank[use][i-1]=b;
				swap++;
			}
		}
	}
	return swap;
}

void presort () {
// called in main to make the initial rankings for scene and world
// view ranks stored as two halves in world->rank
//int 	*data  = Cell::param;
Cell	*scene = Cell::scene;
char	xyz[4] = {'X','Y','Z'};
int	i, j, n = Cell::total;
int	resort = 10;
cellPairs temp[n];
//	data = datain; world = worldin; scene = scenein;
//	sizes = data+N*3; bumps = data+N*4; links = data+N*5; chain = data+N*6;
//	kicks = data+N*7; keeps = data+N*8; repel = data+N*9; bonds = data+N*10;
	printf("Pre-sorting XYZ lists\n"); 
	// pre-sort X
	for (i=0; i<n; i++) { temp[i].a = scene->rank[0][i]; temp[i].d = scene->rank[0][i]->xyz.x; }
	qsort(temp,n,sizeof(cellPairs),sortPairs);
	for (i=0; i<n; i++) { scene->rank[0][i] = scene->rank[0][n+n-i-1] = temp[i].a; }
	// pre-sort Y
	for (i=0; i<n; i++) { temp[i].a = scene->rank[1][i]; temp[i].d = scene->rank[1][i]->xyz.y; }
	qsort(temp,n,sizeof(cellPairs),sortPairs);
	for (i=0; i<n; i++) { scene->rank[1][i] = scene->rank[1][n+n-i-1] = temp[i].a; }
	// pre-sort Z
	for (i=0; i<n; i++) { temp[i].a = scene->rank[2][i]; temp[i].d = scene->rank[2][i]->xyz.z; }
	qsort(temp,n,sizeof(cellPairs),sortPairs);
	for (i=0; i<n; i++) { scene->rank[2][i] = scene->rank[2][n+n-i-1] = temp[i].a; }
	printf("Rough sort done\n"); 
	// re-sort XYZ with added radius shifts for head and tail views
	for (i=0; i<3; i++) {
		printf("Sorting %c\n", xyz[i]);
		for (j=0; j<resort; j++) {
			if (!sortXYZview(i, 1)) break; // sort for head view
		}
		for (j=0; j<resort; j++) {
			if (!sortXYZview(i,-1)) break; // sort for tail view
		}
/*
		// sort local lists
		for (j=0; j<resort; j++) {
			swaps = 0;
			sortXYZcell(world,i);
			if (!swaps) break;
		}
*/
	}
//	data[5] = 0;	// 0 = flag that presort is done
	printf("XYZ lists sorted\n"); 
}
/*
sortXYZrank ( Cells *cell, int use )
{
int	i;
float	ga, gb;
	for (i=1; i<cell->kids; i++)
	{ Cells *a = cell->rank[i-1][use], *b = cell->rank[i][use];
		if (a==0 || b==0) return;
		if (use==0) { ga = a->xyz.x; gb = b->xyz.x; }
		if (use==1) { ga = a->xyz.y; gb = b->xyz.y; }
		if (use==2) { ga = a->xyz.z; gb = b->xyz.z; }
		if (ga > gb) {
			cell->rank[i][use] = a;
			cell->rank[i-1][use]=b;
			a->ranks[use] = i;
			b->ranks[use] = i-1;
			swaps++;
		}
	}
}

sortXYZcell ( Cells *cell, int use )
{
int	i, n, id;
	id = cell->id;
	n = cell->kids;
	if (cell==0) return;
	if (cell->empty) return;
	if (cell->busy < CALM) cell->busy = 0; else cell->busy -= CALM;
	if (n==0) return;
	sortXYZrank(cell,use);
	for (i=0; i<n; i++) sortXYZcell(cell->child[i],use);
	return;
}
*/
