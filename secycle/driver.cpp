#include "sims/util.hpp"
#include "sims/geom.hpp"
#include "sims/cell.hpp"
#include "sims/data.hpp"

void setter () { // called once before driver
	Cell::world->done = 1; // go
}

void driver () { // called on own thread
Cell	*world = Cell::world;
Data	*model = Data::model+world->model;
        if (Cell::world->done == 0) return; // wait
	//Pt(IN driver) Pi(Data::frame) NL
	if (Data::frame > 100) {
		putpdb("model.pdb");
		exit(1);
	}
}

void helper () { // called on same thread as looker
        if (Cell::world->done == 0) return; // wait
	Pt(IN helper) Pi(Data::frame) NL
	sleep(999999);
}

void minder () { // called on same thread as fixers
        if (Cell::world->done == 0) return; // wait
	Pt(IN minder) Pi(Data::frame) NL
	sleep(999999);
}

extern int lastlook;

void looker () { // called on same thread as helper
        if (Cell::world->done == 0) return; // wait
	Pt(IN looker) Pi(Data::frame) NL
	sleep(999999);
}
