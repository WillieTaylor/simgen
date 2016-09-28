#include "util.hpp"
#include "geom.hpp"
#include "cell.hpp"
#include "data.hpp"

void driver ()
{
Cell	*world = Cell::world;
Data	*model = Data::model+world->model;
if (Data::frame < 25000) return;
putpdb("final.pdb");
exit(1);
	if (Data::frame > 5000) {
		putpdb(world);
		exit(1);
	}
	if (Data::frame > 300) Data::tinker = -10;
	if (Data::frame > 100) return;
	DO(i,Cell::total) { Cell *c=Cell::uid2cell[i];
		if (c->level==3) c->xyz.x -= 0.005*(i-Cell::total/2);
	}
	DO(i,Cell::total) { Cell *c=Cell::uid2cell[i];
		if (c->level==2) { Vec x = (c->endC - c->endN)*0.5;
			c->xyz = c->starts->xyz & c->finish->xyz;
			c->endN = c->xyz-x; c->endC = c->xyz+x;
		}
	} 
}
