#include "sims/util.hpp"
#include "sims/geom.hpp"
#include "sims/cell.hpp"
#include "sims/data.hpp"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

void fixRNAstem ( Cell*, float );

static float dmat[300][300];

void setter () { // called once before driver
	FOR(i,Cell::total) { Cell *ci = Cell::uid2cell[i]; int id = ci->atom;
		if (ci->level < depth) continue;
		FOR(j,Cell::total) { Cell *cj = Cell::uid2cell[j]; int jd = cj->atom;
			if (cj->level < depth) continue;
			dmat[id][jd] = dmat[jd][id] = ci->xyz|cj->xyz;
		}
	}
	Cell::world->done = 1; // go
	
}

void minder () { // called on same thread as fixers
	timeout();
        if (Cell::world->done == 0) return; // wait
}

void helper () { // called on same thread as looker
float	d, bump = 1.5, bum = bump*0.5;
int	n = 0;
	timeout();
        if (Cell::world->done == 0) return; // wait
	// if (Data::frame <100) return; // for expansion
	FOR(i,Cell::total) { Cell *ci = Cell::uid2cell[i]; int id = ci->atom;
		if (ci->level < depth) continue;
		FOR(j,Cell::total) { Cell *cj = Cell::uid2cell[j]; int jd = cj->atom;
			if (cj->level < depth) continue;
			if (id==jd) continue;
			if (dmat[id][jd]<3.0 && abs(id-jd)<10) part2cells(ci,cj,dmat[id][jd],0.1);
			if (id==jd+1 || id==jd-1) continue;
			d = ci->xyz|cj->xyz;
			if (d < bum  ) part2cells(ci,cj,bump,-1.0);
			if (d < bump ) part2cells(ci,cj,bump,-0.5);
		}
		if (ci->sis->level==depth && ci->bro->level==depth) {
			d = ci->sis->xyz | ci->bro->xyz;
			if ( d > 2.0 ) part2cells(ci->bro,ci->sis,2.0,0.5);
		}
	}
}

extern int lastlook;

void looker () { // called on same thread as helper
	//timeout();
        if (Cell::world->done == 0) return; // wait
}


void driver () { // called on its own thread
Cell	*world = Cell::world;
int	frame = Data::frame;
Cell	*shell = world->child[0];
Data    *model = Data::model+shell->model;
float   zone = model->sizes[shell->level]*0.5; // radius of shell
int	depth = Data::depth;
int	run = 500;
float	luck = 2.0; // a low value makes links last longer
float	scale = 6.5;
float	d, longest = -1.0;
Cell	*at, *to;
int	it = -1, left = 0;
	timeout();
        if (Cell::world->done == 0) return; // wait
	if (frame >= run) {
		dumpall(scale);
		sortpdb(scale);
		putpdb(world,scale);
		FOR(i,5) { Pi(world->idata[i]) NL }
		exit(1);
	}
/* uncomment for expansion
	shell->xyz.x = shell->xyz.y = shell->xyz.z = 0.0;
	if (frame < 100) { // expand stems
		FOR(i,shell->kids) { Cell *si = shell->child[i]; Vec s = si->xyz;
			if ((si->xyz|shell->xyz) < zone) si->move(s*0.01);
		}
		return;
	}
	if (frame < 200) { // contract
		FOR(i,shell->kids) { Cell *si = shell->child[i]; Vec s = si->xyz;
			if ((si->xyz|shell->xyz) < zone) si->move(-s*0.01);
		}
		return;
	}
*/
//	model->kicks[depth-2] *= 0.99;
	FOR(i,Cell::total) { Cell *ci = Cell::uid2cell[i], *si, *sj;
		if (ci->level < depth) continue;
		FOR(j,ci->nlinks) { Cell *cj = ci->link[j].to;
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
	if (longest < 5.0) return;
	Pr(longest) Pi(at->resn) Pi(to->resn) NL
	// link->next holds the strength (1000..1) +9000 (to stop links breaking when apart)
	d = 1.0 - 0.001*(float)(at->link[it].next - 9000);
	if (at->parent->parent == to->parent->parent) luck *= 10.0; // same stem (inc. chance to break)
	Pi(at->link[it].next-9000) Pr(d) Pr(luck) NL
	if (left==1) { Pt(last link at) Pi(frame) NL } // Pr(d) Pi(at->resn) Pi(at->link[it].to->resn) NL }
	if (drand48() < d*luck+0.01) at->link[it].to = 0;
}
