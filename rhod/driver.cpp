#include "sims/util.hpp"
#include "sims/geom.hpp"
#include "sims/cell.hpp"
#include "sims/data.hpp"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

void pairsout ();
void readlinks ( char*, int );
void tinker ( float, int, int, int );

void driver ()
{
int	frame = Data::frame;
Cell	*world = Cell::world;
Cell	*shell = world->child[0];
Data	*model = Data::model+world->model;
int	depth = Data::depth;
float	rate = 0.0, stop = 993.5; // no shrink
int	run = 1000;
//	WAS float	rate = 1.0, stop = 3.5;
	tinker(0.5,depth,depth,0); // refine just sec.str geom (exc. loops)
	if (frame%10 == 0) { char file[22] = "links pdbid link.dat";
		//putpdb("temp.pdb",world,-Data::scaleout);
		if (model->sizes[1] < 8.0) {
			if (model->sizes[1] < 6.0) {
				if (model->sizes[1]>stop) model->sizes[1] -= 0.02*rate; // 4..5
			} else {
				model->sizes[1] -= 0.05*rate; // 6..8
			}
		} else {
			model->sizes[1] -= 0.10*rate;	 // 8+
		}
		if (frame>=run && model->sizes[1]<stop) {
			DO(i,100) tinker(0.5,depth,depth,0); // refine just sec.str geom (exc. loops)
			pdbout("temp.pdb",world,Data::scaleout);
			pairsout();
			exit(1);
		}
		if (frame > 50) { int links = frame/10;
			Pi(links) NL
			readlinks(file,links); // read increasingly more links
		}
		/*
		DO(i,world->child[0]->kids)
		{ Cell  *sse = world->child[0]->child[i];
		  float s = 0.99;
			sse->endN *= s;
			sse->xyz  *= s;
			sse->endC *= s;
		}
		*/
		/*
		if (frame>500 && frame%100==0) {
			tinker(0.02,depth,depth-1,0); // refine just sec.str geom (exc. loops)
		}
		*/
	}
}

void pairsout () {
// print pairs with separation and depths
// awk '{if($1<$2-10 && $3<2){d=$3+1; s=d*d*(1+$4)*(1+$5); print $1,$2," 0 8 ",s,$3}}' dists.dat | sort -n -k5 > dists.top
int	slot, snap;
int	depth = Data::depth;
FILE	*out = fopen("dists.dat","w");
float	d;
	DO1(i,Cell::total) { Cell *ci = Cell::uid2cell[i];
		if (ci->level != depth) continue;
		DO1(j,Cell::total) { Cell *cj = Cell::uid2cell[j];
			if (cj->level != depth) continue;
			if (ci->resn > cj->resn-12) continue;
			slot = getLink(ci,cj);
			if (slot==0) continue;
			if (slot > 0) {
				snap = ci->link[slot-1].next;
			} else {
				snap = ci->link[1-slot].next;
			}
			d = ci->xyz|cj->xyz;
			fprintf(out,"%d %d %d %f %f %f\n", ci->resn, cj->resn, snap, d, ci->xyz.len(), cj->xyz.len());
		}
	}
}
