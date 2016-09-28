#include "util.hpp"
#include "geom.hpp"
#include "cell.hpp"
#include "data.hpp"

extern int lastlook;

void looker ()
{ // SCALE 0.158 4.75  # scalein = .8/3.8, scaleout = 3.8/.8
FILE	*out = fopen("temp.pdb","w");
Cell	*world = Cell::world;
	DO(i,world->kids)
	{ Cell	*ci = world->child[i];
	  Data	*model = Data::model+ci->model;
	  float	s = Data::scaleout;
		//if (model->moltype==0 && model->subtype==1) s = 3.8/Data::bondCA;
        	putpdb(out,ci,s);
	}
}
