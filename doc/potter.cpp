#include "util.hpp"
#include "geom.hpp"
#include "cell.hpp"
#include "data.hpp"

int computeForces ( Vec*, Vec*, int*, int, int );

void potter ( int level, int cycles ) {
// Run MD refinement at given <level> for specified number of <cycles>
Vec	*pos, *vel, *acc;
int	*cel;
double	dt = 0.001;	// MD time step
float	damp = 0.99;	// velocity damping factor
int	in, n = 0;
int	all = Cell::total;
	pos = new Vec[all];
	vel = new Vec[all];
	acc = new Vec[all];
	cel = new int[all];
	DO1(i,all) { Cell *ci = Cell::uid2cell[i];
		if (ci->level != level) continue;
		pos[n] = ci->xyz;
		vel[n].zero();
		cel[n] = i;
		ci->resn = n;
		n++;
	}
	DO(i,cycles) {
		DO(j,n) {
			vel[j] *= damp;
			pos[j] += vel[j]*dt + acc[j]*0.5*dt*dt;
			vel[j] += acc[j]*0.5*dt;
		}
		in = computeForces(pos,acc,cel,n,level);
		DO(j,n) vel[j] += acc[j]*0.5*dt;
	}
	n = 0;
	DO(i,Cell::total) { Cell *ci = Cell::uid2cell[i];
		if (ci->level != level) continue;
		ci->xyz = (ci->xyz + pos[n++])*0.5;
	}
	Pt(potter) Pi(level) Pi(cycles) Pi(in) NL
}

int computeForces ( Vec *pos, Vec *acc, int *cel, int n, int level ) {
// Gaussian potential with ten times repulsion over attraction (min = -1 at 1)
// plot [0:3] a=14.7378, b=0.801779, c=1/(b*b)
// 10.0*exp(-x*x/(b*b))-a*x*exp(-x*x/(b*b))
// 2*a*c*x*x*exp(-c*x*x)-20*c*x*exp(-c*x*x)-a*exp(-c*x*x)
// force = -(1/r)(dV/dr)
Vec	r;
float	si,sj, ti,tj, ui,uj, rij;
float	a = 14.7378, b =  0.801779;
double	c, d, d2, e, f, p, dp;
int	j, in = 0;
	c = 1.0/(b*b);
	DO(i,n) acc[i].zero();
	for (int i=0; i < n; i++) // for each cell at given level
	{ Cell *ci = Cell::uid2cell[cel[i]], *pai = ci->parent, *cj;
	  Data *pi = Data::model+ci->model, *pj;
		si = pi->sizes[level];
		ti = pi->bonds[level];
		ui = pi->bumps[level];
		if (ci->sis) { // bonded to last
			cj = ci->sis;
			pj = Data::model+cj->model;
			sj = pj->sizes[level];
			tj = pj->bonds[level];
			rij = 0.5*(si+sj+ti+tj);	// bond distance
			r = (ci->xyz-cj->xyz)/rij;	// scale bond distance to 1
			d2 = r.sqr(); d = sqrt(d2);
			p  = 10.0*d2;		// bond deviation = 10d^2
			dp = 10.0*(d-1.0);	// derivitive of quadratic potential
			f = -dp/d;	// divide by d^2 as multiplied by r.xyz (=d) below 
			r *= f;
			j = cj->resn;
			acc[i] += r; acc[j] -= r;
			in++;
		}
		if (ci->bro) { // bonded to last
			cj = ci->bro;
			pj = Data::model+ci->model;
			sj = pj->sizes[level];
			tj = pj->bonds[level];
			rij = 0.5*(si+sj+ti+tj);	// bond distance
			r = (ci->xyz-cj->xyz)/rij;	// scale bond distance to 1
			d2 = r.sqr(); d = sqrt(d2);
			p  = 10.0*d2;		// bond deviation = 10d^2
			dp = 10.0*(d-1.0);	// derivitive of quadratic potential
			f = -dp/d;	// divide by d^2 as multiplied by r.xyz (=d) below 
			r *= f;
			j = cj->resn;
			acc[i] += r; acc[j] -= r;
			in++;
		}
		if (ci->hit) { // bumping
			cj = ci->hit;
			pj = Data::model+cj->model;
			uj = pj->bumps[level];
			rij = 0.5*(ui+uj);		// bump distance
			r = (ci->xyz-cj->xyz)/rij;	// scale bump distance to 1
			d2 = r.sqr();
			if (d2 < 5.0) { // Van der Waal's
				d = sqrt(d2);
				e = exp(-d2*c);
         			p  = e*(10.0-a*d);		// Gaussian based potential (10:1)
				dp = 2.0*a*c*d2*e - 20.0*c*d*e - a*e;	// derivitive of potential
				f = -dp/d;	// divide by d^2 as multiplied by r.xyz (=d) below 
				r *= f;
				j = cj->resn;
				acc[i] += r; acc[j] -= r;
				in++;
			}
		}
	}
	return in;
}
