#include "sims/util.hpp"
#include "sims/geom.hpp"
#include "sims/cell.hpp"
#include "sims/data.hpp"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

void fixRNAstem ( Cell*, float );

void driver ()
{
int	frame = Data::frame;
Cell	*world = Cell::world;
Cell	*shell = world->child[0];
Data    *model = Data::model+shell->model;
float   zone = model->sizes[shell->level]*0.5; // radius of shell
int	depth = Data::depth;
int	run = 300;
float	luck = 2.0; // a low value makes links last longer
float	scale = 3.5;
float	d, longest = -1.0;
Cell	*at, *to;
int	it = -1, left = 0;

if (frame >= run) { sortpdb(scale); exit(1); }

	model->kicks[depth-2] *= 0.99;
	DO(i,Cell::total) { Cell *ci = Cell::uid2cell[i], *si, *sj;
		if (ci->level < depth) continue;
		DO(j,ci->nlinks) { Cell *cj = ci->link[j].to;
			if (cj == 0) continue;
			left++;
			si = ci->parent->parent;
			sj = cj->parent->parent;
			if ((si->xyz|si->parent->xyz) > zone) continue;
			if ((sj->xyz|sj->parent->xyz) > zone) continue;
			// cull only when stems are both inside the shell diam
			d = ci->xyz | cj->xyz;
			if (d < longest) continue;
			it = j; at = ci; to = cj;
			longest = d;
		}
	}
	if (it<0) return;
	//Pr(longest) Pi(at->resn) Pi(to->resn) NL
	// link->next holds the strength (1000..1) +9000 (to stop links breaking when apart)
	d = 1.0 - 0.001*(float)(at->link[it].next - 9000);
	if (at->parent->parent == to->parent->parent) luck *= 10.0; // same stem (inc. chance to break)
	//Pi(at->link[it].next-9000) Pr(d) Pr(luck) NL
	if (left==1) { Pt(last link at) Pi(frame) NL } // Pr(d) Pi(at->resn) Pi(at->link[it].to->resn) NL }
	if (drand48() < d*luck+0.01) at->link[it].to = 0;
	if (left==1 || (frame >= run)) {
		fixRNAstem(Cell::world,1.0);
		sortpdb(scale);
		exit(1);
	}
}
