#include "sims/util.hpp"
#include "sims/geom.hpp"
#include "sims/cell.hpp"
#include "sims/data.hpp"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

void driver ()
{
int	frame = Data::frame;
Cell	*world = Cell::world;
Cell	*shell = world->child[0];
Data	*model = Data::model+world->model;
int	depth = Data::depth;
int	total = Cell::total;
int	run = 2000;
float	scale = 1.0;
Vec	*a = new Vec[total/2];
Vec	*b = new Vec[total/2];
Cell	*p0 = world->child[0];
Cell	*p1 = world->child[1];
Vec	in = Vec(0.01,0,0);
int	m = 0, n = 0;
int	aa,bb,ab;
float	d, bump = 0.5;
        if (frame < 100) return;
	if (frame == 100) pdbout("start.out");
	if (frame%10==0) { // count bumps
		FOR(i,total) { Cell *c = Cell::uid2cell[i], *p; int k;
			if (c->level < depth) continue;
			p = c->parent->parent;
			if (depth==3) k = p->id;
			if (depth==4) k = p->parent->id;
			if (depth==5) k = p->parent->parent->id;
			if (depth==6) k = p->parent->parent->parent->id;
			if (depth==7) k = p->parent->parent->parent->parent->id;
			if (k==0) a[m++] = c->xyz;
			if (k==1) b[n++] = c->xyz;
		}
		aa = 0;
		FOR(i,m) {
			FOR(j,m) {
				if (i-j < 2) continue;
				d = a[i]|a[j];
				if (d>bump) continue;
				aa++;
			}
		}
		bb = 0;
		FOR(i,n) {
			FOR(j,n) {
				if (i-j < 2) continue;
				d = b[i]|b[j];
				if (d>bump) continue;
				bb++;
			}
		}
		ab = 0;
		FOR(i,m) {
			FOR(j,n) {
				d = a[i]|b[j];
				if (d>bump) continue;
				ab++;
			}
		}
		Pi(frame) Pi(aa) Pi(bb) Pi(ab) NL
	}
	p0->move(-in); p1->move(in);
	if (p0->xyz.x<0 && p1->xyz.x>0) {
	//if (p0->xyz.x<8 && p1->xyz.x>8) {
	//if (frame > 1000) {
		pdbout("final.out");
		pdbout("protA.out",p0);
		pdbout("protB.out",p1);
		exit(1);
	}
}
