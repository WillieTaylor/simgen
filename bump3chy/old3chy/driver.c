#include "main/common.h"

void driver ( int *datain, Cells *worldin )
{
int	i, wait = 1000, frame, run = 100000;
float	speed = 0.50;
	run = (int)(0.01/speed*(float)run);
	data = datain; world = worldin;
	uid2cell = world->link;
	frame = data[0];
	total = data[1];
	depth = data[2];
	natoms = data[8];
	if (frame < wait) return;
	if (frame==wait) {
		Pt(GO) NL
	}
	if (frame%run > run/2) { 
		world->child[0]->xyz.x += speed;
		world->child[1]->xyz.x -= speed;
		world->child[0]->xyz.y *= 0.99;
	} else {
		world->child[0]->xyz.x -= speed;
		world->child[1]->xyz.x += speed;
		world->child[1]->xyz.y *= 0.99;
	}
}
