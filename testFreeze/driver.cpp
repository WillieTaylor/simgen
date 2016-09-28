#include "sims/util.hpp"
#include "sims/geom.hpp"
#include "sims/cell.hpp"
#include "sims/data.hpp"

void driver ()
{
Cell	*world = Cell::world;
Data	*model = Data::model+world->model;
	world->child[0]->child[0]->xyz.x -= 0.1;
}
