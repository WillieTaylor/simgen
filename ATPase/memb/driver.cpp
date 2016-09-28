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
Cell	*A = world->child[0];
Cell	*C = world->child[1];
Data	*model = Data::model+world->model;
int	depth = Data::depth;
int	total = Cell::total;
int	run = 2000;
float	spin = 0.5;
	C->spin(C->endN,C->endC,spin*(randf()-0.5));
}
