#include "util.hpp"

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include "geom.hpp"
#include "cell.hpp"
#include "data.hpp"

#define N 10	// levels of granularity for each object

typedef struct {
	int	type, level, sort;
	GLuint	*body;	// points to a series of N bodies at different levels of detail
	Vec	endN, endC, axes;
	int	flip;
	Mat	rot;
} Types;
static Types ***types;	// holds an array of sorts for each types for each model

static GLint	fogMode;
static GLfloat  dout = 100.0; // eye position on Z
static GLfloat  back = 500.0; // backplane cutoff 
static GLfloat  tx, ty, tz;
static GLfloat  rx, ry, rz;
static Mat	R, View;
static Vec	eye, focus;

static float trX, trY, trZ;
static float rtX, rtY, rtZ;

static float maxthick = 5;
static GLuint sticks, tube[6], mesh[6], bond[6];
static GLfloat ***colour;
static float see[3];

void sceneCell ( Cell *cell, Cell *scene )
{
int     n = cell->kids;
        if (n==0) return;
        FOR(i,n) {
                FOR(j,3) scene->rank[j][total] = scene->rank[j][total+Cell::total] = cell->child[i];
		total++;
        }
        FOR(i,n) sceneCell(cell->child[i], scene);
}

Cell* setScene ( Cell *world )
{
int	n = Cell::total * 2;
Cell	*scene;
	scene = new Cell;
        scene->rank = new Cell**[3]; TEST(scene->rank)
        FOR(j,3) { scene->rank[j] = new Cell*[n]; TEST(scene->rank[j]) }
	total = 0;
        sceneCell(world,scene);     // sets-up front and rear views
	return scene;
}

void hemiSphere(double rad, int gran, float flip ) {
GLUquadricObj* quad = gluNewQuadric(); 
float	g = (float)gran, r0, r1, d, dr = 2.0*rad/g, t, dt = halfPI/g;
	r0 = rad; t = 0.0;
	FOR(i,gran) {
		t += dt;
		d = dr*cos(t);
		r1 = rad*cos(t);
		if (flip>0.0) {
			gluCylinder( quad, r0, r1, d, gran, 1 );
			glTranslatef( 0.0, 0.0, d*flip);
		} else {
			glTranslatef( 0.0, 0.0, d*flip);
			gluCylinder( quad, r1, r0, d, gran, 1 );
		}
		r0 = r1;
	}
}

void makeBody ( Types*, int, int, int );

static float conect = 0.3;
static int   meshed = 8;

static int last=0, hold=0;

int	detail[N+1] = {0,30,28,26,24,22,20,18,16,14,12}; // starting detail at each level

GLuint makeTubes ( float thick, float length, float scale, int cap, int grain, int solid )
{
// cap: 0 = open, 1 = one cap, 2 = two caps. NB Z not scaled
GLuint	cylinder;
GLUquadricObj* quad = gluNewQuadric(); 
GLuint	gran = (GLuint)grain;
GLdouble rad = (GLdouble)thick*0.5;
GLdouble len = (GLdouble)length;
float	flip = 1.0;
	if (cap < 0) { cap = -cap; flip = -1.0; }
	cylinder = glGenLists(1);
	glNewList(cylinder, GL_COMPILE);
		glPushMatrix();
			glTranslatef( 0.0, 0.0, -0.5*len);
			glScalef(scale, scale, 1.0);
			if (solid<0) {
				glPolygonMode ( GL_FRONT_AND_BACK, GL_LINE );
			} else {
				glPolygonMode ( GL_FRONT_AND_BACK, GL_FILL );
			}
			gluCylinder( quad, rad, rad, len, gran, gran );
		glPopMatrix();
		if (cap > 0) { // cap one end
			glPushMatrix();
				glTranslatef( 0.0, 0.0, 0.5*len);
				glScalef(scale, scale, 1.0);
				if (solid<0) {
					glPolygonMode ( GL_FRONT_AND_BACK, GL_LINE );
				} else {
					glPolygonMode ( GL_FRONT_AND_BACK, GL_FILL );
				}
				hemiSphere(rad, gran, 1.0);
			glPopMatrix();
		}
		if (cap > 1) { // cap other end
			glPushMatrix();
				glTranslatef( 0.0, 0.0, -0.5*length);
				glScalef(scale, scale, 1.0);
				if (solid<0) {
					glPolygonMode ( GL_FRONT_AND_BACK, GL_LINE );
				} else {
					glPolygonMode ( GL_FRONT_AND_BACK, GL_FILL );
				}
				hemiSphere(rad, gran,-1.0);
			glPopMatrix();
		}
	glEndList();
	return cylinder;
}

void makeTypes ()
{
/*
 * each model can have a different set of body shapes.
 * each level has a shape type (sphere/tube/ellipsoid) specfied in Data.
 * different levels can have a different size (eg; big/small spheres).
 * each type can have different sorts (eg; diameter:length) set by sub[]
 * each sort has a series of bodies (body[]) at different resolutions
 * all referred to through types[model][type][sort].body[n]
 *
 * type 
 * 0 = tiny sphere
 * 1 = sphere
 * 2 = tube (different caps)
 * 3 = ellipsoid
 */
int	i;
// sorts for type: 0 1 2 3 4 5 6 7 8 9  (4-9 not used yet)
int	sorts[N] = { 1,1,10,E,0,0,0,0,0,0};
int	level, type, k;
int	nmodels = Data::nmodels;
	// types allocation for nmodels made in startup()
	for (i=0; i<nmodels; i++) {
		sizes = Data::model[i].sizes;
		types[i] = (Types**)malloc(sizeof(Types*)*N); TEST(types[i]) // allocate for sorts
		types[i][0] = 0;	// level 0 = world     not rendered
		for (level=1; level<N; level++) // make types for each level
		{ int	kind = Data::model[i].shape[level];
			k = abs(kind);
			types[i][level] = (Types*)malloc(sizeof(Types)*sorts[k]); TEST(types) // allocate for sorts
			makeBody(types[i][level],level,kind,sorts[k]);
		}
	}
}

void makeBody ( Types *type, int level, int kind, int sorts )
// Make body variations for each level
{
float	length, radius, angle, linkto, len, dx,dy,dz, cosA,sinA, d2r = PI/180.0;
GLuint	body, link;
int	i,j,k, kinds, solid, kdrop = 2;
//float	*sizes = Data::model[0].sizes;
	if (kind > 0) { kinds =  kind; solid =  1; }
	if (kind ==0) { kinds =  kind; solid =  1; } // virtual = solid?
	if (kind < 0) { kinds = -kind; solid = -1; }
	for (i=0; i<sorts; i++) { // same set-up for each type
		type[i].sort = i;
		type[i].type = kind;
		type[i].level = level;
		type[i].body = (GLuint*)malloc(sizeof(GLuint)*N); TEST(type[i].body)
		type[i].axes.x = 1.0;
		type[i].axes.y = 1.0;
		type[i].axes.z = 1.0;
		type[i].endN.zero();
		type[i].endC.zero();
		type[i].rot = Mat();
	}
	// the following scale the ellipsoid Z axis from 2/E to E/2 over E intervals

	switch (kinds) {

		case 0 : // virtual sphere (rendered as tiny sphere)
			for (i=0; i<sorts; i++)	{	// make each variant (just 1)
				k = detail[level];	// starting level of detail (drops by 2 each time)
				for (j=0; j<N; j++) {	// make each level of detail
					type[i].body[j] = glGenLists(1);
					glNewList(type[i].body[j], GL_COMPILE);
						glutSolidSphere(0.1, 6, 6);
					glEndList();
					k -= kdrop; if (k<3) k = 3;
				}
				type[i].flip = 0;
			}
			// endN and endC are set for each cell in getFrame()
		return;

		case 1 : // sphere
			radius = 0.5;
			for (i=0; i<sorts; i++)		// make each variant (just 1)
			{ GLdouble rad = (GLdouble)radius;
				k = detail[level];	// starting level of detail (drops by 2 each time)
				for (j=0; j<N; j++) {	// make each level of detail
					type[i].body[j] = glGenLists(1);
					glNewList(type[i].body[j], GL_COMPILE);
						if (solid > 0) {
							glutSolidSphere(rad, k, k);
						}
						if (solid < 0) {
							glutWireSphere( rad, k, k);
						}
					glEndList();
					k -= kdrop; if (k<3) k = 3;
				}
				type[i].flip = 0;
			}
			// endN and endC are set for each cell in getFrame()
		return;

		case 2 : // tube = cylinder + caps + links
			for (i=0; i<sorts; i++) {	// make each sort with different lengths
				radius = 0.5;
				length = radius*(1.0+(float)i); // rad:len = 1:1 ... 1:<sorts> 
				k = detail[level];	// starting level of detail (drops by 2 each time)
				for (j=0; j<N; j++) {		// make each level of detail
					body = makeTubes(radius*2.0, length, 1.0, 2, k, solid);
					type[i].body[j] = glGenLists(1);
					glNewList(type[i].body[j],  GL_COMPILE);
						glPushMatrix();
							glCallList(body);
						glPopMatrix();
					glEndList();
					k -= kdrop; if (k<4) k = 4;
				}
				dz = length*0.5;
				dz += radius*0.5; // shift ends more into nose
				type[i].endN.z = -dz;
				type[i].endC.z =  dz;
				type[i].flip =  2;	// when sort==2 body will lie along N-1 --> N+1
			}
		return;

		case 3 : // ellipsoid (connects at poles), ratio range: 1=0.1, oblate 10=1, prolate 19=10
			radius = 0.5;
			for (i=1; i<sorts; i++) { // NL runs 1...E-1 (type[0] not set)
				dx = type[i].axes.x = 1.0; 
				dy = type[i].axes.y = 1.0; 
				dz = type[i].axes.z = Data::Eratio[i];
				k = detail[level];	// starting level of detail (drops by 2 each time)
				for (j=0; j<N; j++) {	// make each level of detail
					type[i].body[j] = glGenLists(1);
					glNewList(type[i].body[j], GL_COMPILE);
						glPushMatrix();	// make body
							glScalef(dx, dy, dz);
							if (solid > 0) {
								glutSolidSphere(radius, k, k);
							}
							if (solid < 0) {
								glutWireSphere(radius, k, k);
							}
						glPopMatrix();
/*
						if (chain[level]) {
							glPushMatrix();	// make connector
								glTranslatef( 0.0, 0.0,  dz*radius);
								if (solid > 0) {
									glutSolidSphere(conect, 8, 8);
								}
								if (solid < 0) {
									glutWireSphere( conect, 8, 8);
								}
							glPopMatrix();
							glPushMatrix();	// make connector
								glTranslatef( 0.0, 0.0, -dz*radius);
								if (solid > 0) {
									glutSolidSphere(conect, 8, 8);
								}
								if (solid < 0) {
									glutWireSphere( conect, 8, 8);
								}
							glPopMatrix();
						}
*/
					glEndList();
					k -= kdrop; if (k<3) k = 3;
				}
				dz *= radius;
				type[i].endN.x = 0.0; type[i].endN.y = 0.0; type[i].endN.z = -dz;
				type[i].endC.x = 0.0; type[i].endC.y = 0.0; type[i].endC.z =  dz;
				type[i].flip = 2; // Z will align with N-1 --> N+1
			}
		return;
	}
}

static void init ( void )
//  Initialize depth buffer, fog, light source, material property, and lighting model.
{
GLfloat position[4] = { 0.5, 0.5, 3.0, 0.0 };
int	i, j, k;
double	size;
	see[0] = 0.05; see[1] = 0.3; see[2] = 1.0; // 0=cthru, 1=almost, 2=solid
	colour = (GLfloat***)malloc(sizeof(GLfloat**)*3); TEST(colour)
	for (i=0; i<3; i++) { float set = 0.0;
		colour[i] = (GLfloat**)malloc(sizeof(GLfloat*)*N); TEST(colour[i])
		for (j=0; j<N; j++) {
			if (j>2) set = 1.0;
			colour[i][j] = (GLfloat*)malloc(sizeof(GLfloat*)*4); TEST(colour[i][j])
			for (k=0; k<3; k++) colour[i][j][k] = set; 
			colour[i][j][3] = see[i];
		}
		//      = transparency
		//     |   = preset (col)
		//     |  |   = RGB (0,1,2)
		//     |  |  |
		colour[i][0][2] = 1.0; // blue
		colour[i][1][1] = 1.0; // green
		colour[i][2][0] = 1.0; // red
		colour[i][3][0] = 0.0; // cyan
		colour[i][4][1] = 0.0; // magenta
		colour[i][5][2] = 0.0; // yellow
		colour[i][6][0]=colour[i][6][1]=colour[i][6][2]=0.5; // grey
		colour[i][7][0]=colour[i][7][1]=colour[i][7][2]=0.0; // black
	}			// rest white
	// colour[i=123][N][012=rgb,3=see[i]], see = 0=cthru, 1=almost solid, 2=solid
	//glEnable(GL_RESCALE_NORMAL);
	glEnable(GL_NORMALIZE);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable (GL_BLEND);

	glLightfv(GL_LIGHT0, GL_POSITION, position);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	{ GLfloat mat[4];
		mat[0] = 0.0; mat[1] = 0.0; mat[2] = 0.0;		
		glMaterialfv(GL_FRONT, GL_AMBIENT, mat);
		mat[0] = 1.0; mat[1] = 1.0; mat[2] = 1.0;		
		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat);
		mat[0] = 1.0; mat[1] = 1.0; mat[2] = 1.0;
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat);
		glMaterialf (GL_FRONT, GL_SHININESS, 50.0);
	}
	glEnable(GL_FOG);
	{ GLfloat fogColour[4] = {0.5, 0.5, 0.5, 1.0};
	  float	fade = 1.0;
		fogMode = GL_EXP2;
		//fogMode = GL_LINEAR;
		glFogi (GL_FOG_MODE, fogMode);
		glFogfv(GL_FOG_COLOR, fogColour);
		glFogf (GL_FOG_DENSITY, 0.003); //0.005);
		glHint (GL_FOG_HINT, GL_DONT_CARE);
		glFogf (GL_FOG_START, dout+fade);	// only active 
		glFogf (GL_FOG_END, dout+fade+fade);	// with LINEAR
	}
	glClearColor(0.5, 0.5, 0.5, 1.0); // should be same as fog colour

	tx = ty = tz = 0.0;
	rx = ry = rz = 0.0;
	trX = trY = trZ = 0.0;
	rtX = rtY = rtZ = 0.0;

	tube[0] = makeTubes(conect, 1.0, 0.3, 0, meshed, 1);
	tube[1] = makeTubes(conect, 1.0, 0.5, 0, meshed, 1);
	tube[2] = makeTubes(conect, 1.0, 0.7, 0, meshed, 1);
	tube[3] = makeTubes(conect, 1.0, 0.9, 0, meshed, 1);
	tube[4] = makeTubes(conect, 1.0, 1.1, 0, meshed, 1);
	tube[5] = makeTubes(conect, 1.0, 1.3, 0, meshed, 1);

	mesh[0] = makeTubes(conect, 1.0, 0.3, 0, meshed,-1);
	mesh[1] = makeTubes(conect, 1.0, 0.5, 0, meshed,-1);
	mesh[2] = makeTubes(conect, 1.0, 0.7, 0, meshed,-1);
	mesh[3] = makeTubes(conect, 1.0, 0.9, 0, meshed,-1);
	mesh[4] = makeTubes(conect, 1.0, 1.1, 0, meshed,-1);
	mesh[5] = makeTubes(conect, 1.0, 1.3, 0, meshed,-1);

	bond[0] = makeTubes(conect, 1.0, 0.1, 0, meshed, 1);
	bond[1] = makeTubes(conect, 1.0, 0.2, 0, meshed, 1);
	bond[2] = makeTubes(conect, 1.0, 0.3, 0, meshed, 1);
	bond[3] = makeTubes(conect, 1.0, 0.4, 0, meshed, 1);
	bond[4] = makeTubes(conect, 1.0, 0.5, 0, meshed, 1);
	bond[5] = makeTubes(conect, 1.0, 0.6, 0, meshed, 1);
}

void getFrame ( Cell *cell )
{
Vec	x,y,z;
	if (cell->level==Data::depth) return;	// assume atoms are spheres
	x = cell->endC - cell->endN;		// x = tube/ellipsoid axis
	y = cell->xyz;				// 'random' y since axial symmetry
	x.setVec(); y.setVec();			// normalise
	z = (x^y).getNorm(); y = (z^x).getNorm(); // orthogonalise + normalise
	cell->rot = Mat(y,z,x);			// assumes axis is along Z
}

void orientBody ( Cell *cell )
{
int	i, id = cell->id, n = cell->kids;
	if (n==0) return;
	if (cell->empty) return;
	for (i=0; i<n; i++)
	{ Cell *child = cell->child[i];
		if (child->empty) continue;
		getFrame(child);
		orientBody(child);
	}
}

void fillCell (Cell *cell)
{
int     ii, i, j, k, n = cell->kids;
GLfloat col[4] = {0.0, 0.0, 0.0, 0.8};
float   a,A, b,B, c,C, d = 0.03*(float)sizes[1];
float   norm = 0.5, flick = 1.0, sumup = 0.0;
        for (i=0; i<n; i++) { Vec u,v,w;
                if (i) j = i-1; else j = n-1;
		u = cell->child[i]->xyz - cell->xyz;
                v = cell->child[j]->xyz - cell->xyz;
                w = v^u;
                sumup += w.z;
        }
        if (cell->model==0) col[1] = 1.0;
        if (cell->model==1) col[0] = 1.0;
        if (cell->model==2) col[2] = 1.0;
        glMaterialfv(GL_FRONT, GL_DIFFUSE, col);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        for (ii=0; ii<n; ii++) {
                i = ii;
                if (i) j = i-1; else j = n-1;
                if (sumup < 0.0) { k=i; i=j; j=k; }
                glBegin (GL_TRIANGLES);
                        glNormal3f(0.0, 0.0, norm+randf()*flick);
                        glVertex3f(cell->child[j]->xyz.x,cell->child[j]->xyz.y,cell->child[j]->xyz.z);
                        glVertex3f(cell->xyz.x, cell->xyz.y, cell->xyz.z + d);
                        glVertex3f(cell->child[i]->xyz.x,cell->child[i]->xyz.y,cell->child[i]->xyz.z);
                glEnd();
                glBegin (GL_TRIANGLES);
                        glNormal3f(0.0, 0.0, -norm*randf()*flick);
                        glVertex3f(cell->child[j]->xyz.x,cell->child[j]->xyz.y,cell->child[j]->xyz.z);
                        glVertex3f(cell->xyz.x, cell->xyz.y, cell->xyz.z - d);
                        glVertex3f(cell->child[i]->xyz.x,cell->child[i]->xyz.y,cell->child[i]->xyz.z);
                glEnd();
        }
}

void drawLink ( Vec end1, Vec end2, int type, int sort )
{ // type: 1=sphere, 2=tube, 3=ellipsoid. sort: 0 = bond, >0 = link, <0 = mesh
float	angle, len;
Vec	x,y,z;
Vec	mid = end1 & end2;
	len = end1 | end2;
	if (sort > maxthick) sort = maxthick;
	z.x = 0.0; z.y = 0.0; z.z = 1.0;
	y.setVec(end1-end2);
	angle = acos(z*y)*180.0/PI;
	x = z^y;
	glPushMatrix();
		glTranslatef(mid.x, mid.y, mid.z);
		glRotatef(angle, x.x, x.y, x.z);
		glScalef(1.0,1.0,len);
		if (type > 0) glCallList(tube[sort]);
		if (type ==0) glCallList(bond[sort]);
		if (type < 0) glCallList(mesh[sort]);
	glPopMatrix();
}

void drawBonds ( Cell *cell, Cell *bond, int bondlink) {
//	bondlink = cell->bond[bid].link: 1=NN, 2=NC, 3=CN, 4=CC (-ve = crosslink)
int	mod = cell->model;
Data	*param = Data::model+mod;
float	size = param->sizes[cell->level];
int	moltype = param->moltype,
	subtype = param->subtype;
Cell	*parent = cell->parent;
//Bonds	*cbonds = cell->bond;
int	i, j, m, id,
	level = cell->level,
	solid = cell->solid,
	type = abs(cell->type),
	sort = cell->sort,
	bondto[4], nbonds,
	bondtype;
float	hid[3] = { loopTHIC, alphTHIC, betaTHIC },
	sunk = size*(0.5-0.05); // radius less a bit  
Vec	v, n, c,
	endN = cell->endN,
	endC = bond->endC;
int	thick = depth-level;
	if (type==0) return; // don't bond virtual centres
	if (cell->type < 0) thick = 0;// thin links for wireframes
	if (bondlink < -990) return; // don't draw a RELINKed bond
	if (moltype==1 && level<depth && cell->type==2) return; // no nucleic tube bonds
	if (level != bond->level) return;
	if (bondlink==0 && cell->type==1 && bond->type==1) { // sphere--sphere
		v = (bond->xyz - cell->xyz).getVec(sunk);
		n = cell->xyz + v;
		c = bond->xyz - v;
		drawLink(n, c, type, thick);
		return;
	}
	if (bondlink > 0) { // specified pole--pole connection
		if (bondlink==4 || bondlink==3) endN = cell->endC;
		if (bondlink==1 || bondlink==3) endC = bond->endN;
		if (type==1) { // sphere
			v = (bond->xyz - cell->xyz).getVec(sunk);
			n = cell->xyz + v;
			c = bond->xyz - v;
		}
		if (type==2) { // tube
			v = (endC-endN).getVec(sunk);
			n = endN + v*hid[sort];
			c = endC - v*hid[bond->sort];
		}
		if (type==3) { // ellipsoid
			n = endN;
			c = endC;
		}
		drawLink(n, c, type, thick);
		return;
	}
	if (bondlink==0) // link closest poles
	{ float dnn = vdif(cell->endN,bond->endN),
		dnc = vdif(cell->endN,bond->endC),
		dcn = vdif(cell->endC,bond->endN),
		dcc = vdif(cell->endC,bond->endC),
		limit = 1.2*(float)sizes[level];
	  	bondtype = sort4min(dnn,dnc,dcn,dcc);
	  	switch (bondtype) {
			case 0 : // if all equal, default = N--C
				if (dnc > limit) return;
				drawLink(cell->endN, bond->endC, type, thick);
			break;
			case 1 :
				if (dnn > limit) return;
				drawLink(cell->endN, bond->endN, type, thick);
			break;
			case 2 :
				if (dnc > limit) return;
				drawLink(cell->endN, bond->endC, type, thick);
			break;
			case 3 :
				if (dcn > limit) return;
				drawLink(cell->endC, bond->endN, type, thick);
			break;
			case 4 :
				if (dcc > limit) return;
				drawLink(cell->endC, bond->endC, type, thick);
			break;
		}
		return;
	}
	if (bondlink < 0) { // ?
		bondlink = -bondlink;
Px
exit(1);
	}
}
//	model = cell->model;
//	m = model*M*N+N;
//	moltype = data[m];
//	links = data+m+N*5;
//	chain = data+m+N*6;
//        nbonds = abs(chain[level]);
//	if (moltype==2) { // CHEMical bonds
//		if (level==depth && parent->solid<=0) {
//			for (i=0; i<nbonds; i++) {
//				if (cbonds[i].to==0) continue; // don't break as cyclic bonds are at the end
//				if (cbonds[i].to->level != level) continue;
//				if (cell->uid > cbonds[i].to->uid) continue;
//				if (cbonds[i].next<0) type = -type; // -ve = mesh to mark cyclic pseudo-bonds
//				drawLink(cell->xyz, cbonds[i].to->xyz, type, cbonds[i].type);
//			}
//		}
//	} else { int bondtype;
//		if (cell->bond && cell->bond[0].type ) {
//			bondtype = depth-level;
//			for (i=0; i<5; i++) { Bonds *hinge = cell->bond+i;
//				if (hinge->to == 0) break;
//				if (hinge->next==1)  drawLink(cell->endN, hinge->to->endN, type, bondtype);
//				if (hinge->next==2)  drawLink(cell->endN, hinge->to->endC, type, bondtype);
//				if (hinge->next==3)  drawLink(cell->endC, hinge->to->endN, type, bondtype);
//				if (hinge->next==4)  drawLink(cell->endC, hinge->to->endC, type, bondtype);
//			}
//		}
//		if (cell->sis && cell->sis->level==cell->level) { Cells *bond = cell->sis; int inseq;
//			bondtype = 1;
//			if (level==depth) {
//				if (abs(type)==1) { // atom spheres (-ve type = mesh)
//					drawLink(cell->xyz, bond->xyz, type, bondtype);
//				} else { // for other shapes
//					drawLink(cell->endN, bond->endC, type, bondtype);
//				}
//				return;
//			}
//			if (moltype==0) { // protein
//				if (type==1) inseq = 1;		// spheres
//				if (type==3) inseq = -1;	// ellipsoids
//				if (chain[level]<0) inseq = -1;	// cyclic
//				if (level==depth-1) inseq = 1;	// SSEs
//			} else {	// RNA
//				if (level<depth) inseq = -1;	// SSEs or domains
//			}
//			if (inseq==1) { // join ends in sequence
//				drawLink(cell->endN, bond->endC, type, bondtype);
//				return;
//			} else	     // join ends by proximity
//			{ float dnn = vdif(cell->endN,bond->endN),
//				dnc = vdif(cell->endN,bond->endC),
//				dcn = vdif(cell->endC,bond->endN),
//				dcc = vdif(cell->endC,bond->endC),
//				limit = 0.2*(float)sizes[level];
//				//if (cell->prox.x < NOISE) limit = -9.9; // unrefined bond (thin)
//				if (cell->type < 0) limit = -9.9;	// thin links for wireframes
//			  	switch (fsort4min(dnn,dnc,dcn,dcc)) {
//					case 0 :
//						if (dnc > limit) bondtype = 0;
//						drawLink(cell->endN, bond->endC, type, bondtype);
//					break;
//					case 1 :
//						if (dnn > limit) bondtype = 0;
//						drawLink(cell->endN, bond->endN, type, bondtype);
//					break;
//					case 2 :
//						if (dnc > limit) bondtype = 0;
//						drawLink(cell->endN, bond->endC, type, bondtype);
//					break;
//					case 3 :
//						if (dcn > limit) bondtype = 0;
//						drawLink(cell->endC, bond->endN, type, bondtype);
//					break;
//					case 4 :
//						if (dcc > limit) bondtype = 0;
//						drawLink(cell->endC, bond->endC, type, bondtype);
//					break;
//				}
//			}
//		}
//	}
//}

void drawLinks ( Cell *cell ) {
Cell	*parent = cell->parent;
//Bonds	*cbonds = cell->bond;
int	i, j, m, id,
	level = cell->level,
	solid = cell->solid,
	type = cell->type,
	sort = cell->sort,
	model = cell->model,
	bondto[4], nbonds;
Data	*param = Data::model+model;
	moltype = param->moltype;
	links = param->links;
	chain = param->chain;
	// draw links in a chain (drawLink(...,0) = thin)
	if (cell->link == 0) return;
	for (i=0; i<links[level]; i++)
	{ Cell *linki = cell->link[i].to;
	  int	thick = cell->link[i].type;
		if (linki == 0) continue;
		if (linki->empty) {
			//cell->link[i] = 0;
			continue;
		}
		thick = depth-level;
		drawLink(cell->xyz,linki->xyz,0,thick);
	}
}

float setCthru ( Cell *cell ) {
// calculate opacity of the cell from its distance/(angle) from the viewer (eye)
Cell	*parent = cell->parent;
float	profile[5] = { 5.0, 5.0, 5.0, 5.0, 5.0 };
float	aa, bb, cc, cosm, backcut = back*back;
int	solid = cell->solid,
	level = cell->level,
	type  = cell->type;
Vec	cent;
float	cthru, range = 0.00001,
	hide = (float)Data::hidden;
	if (hide > 0) range /= hide;	// becomes solid farther
	if (hide < 0) range *= -hide;	// becomes solid closer
        if (solid < -99) cell->solid = -99;
        if (solid >  99) cell->solid =  99;
        if (solid < parent->solid) cell->solid = parent->solid;
	if (cell->heat > 0.5) cell->solid = -999;     // always see into an active object
        if (type <= 0) {
                cell->solid = -999;     // always see into a wireframe or virtual cell
                type = -type;
        }
        bb = eye||cell->xyz;		// NB bb = b**2
        if (bb > backcut) {		// out of sight
                cell->solid++;
                return -1.0; // -ve = flag to skip rendering
        }
	cent = focus;
        aa = cent||cell->xyz;
        bb =  eye||cell->xyz;
        cc =  eye||cent;
        cell->far = 0.5*(cell->far + bb);               // smoothly shift from last 
        cthru = see[0]+1.0-exp(-cell->far*range);	// become more opaque with distance
	// angle condition
        cosm = 0.5*(bb+cc-aa)/(sqrt(bb)*sqrt(cc));      // cos=1 is dead ahead
        cc = 1.0 - cosm;                                // ahead = 0, behind = 2
        cc = 2.0 - exp(-cc*cc*5.0);                     // ahead = 1 rising to 2
        if (cosm < 0.7 ) {                              // exclude outside +/-45 deg. cone
                cell->solid++;
                return -1.0;
        }
        if (level==depth) cthru = 1.0;
        if (level==depth-1) cthru *= profile[0];
        if (level==depth-2) cthru *= profile[1];
        if (level==depth-3) cthru *= profile[2];
        if (level==depth-4) cthru *= profile[3];
        if (level==depth-5) cthru *= profile[4];
        if (cthru > 0.9) {	// if transparency low, move towards solid
                cell->solid++;
        } else {
                if (cell->live && parent->solid<0) cell->solid--;
        }
	return cthru;
}

static void drawCell (Cell *cell)
{
Cell	*parent = cell->parent,
	*hit = cell->hit;
int	i, m, id, in, gran, lost,
	model = cell->model,
	level = cell->level,
	solid = cell->solid,
	type = cell->type,
	sort = cell->sort,
	rat, col, colev, cols[10] = { 6, 1, 2, 3, 4, 5, 0, 1, 2, 3 };
	// 0=blue, 1=green, 2=red, 3=cyan, 4=magenta, 5=yellow, 6=grey, 7=black, 8=white
float	heat = cell->heat;
Vec	link, last;
Data	*param = Data::model+model;
int	nbonds = abs(param->chain[level]);
float	r, s, t, size = param->sizes[level];
int	moltype = param->moltype,
	subtype = param->subtype,
	colours = param->colours;
GLfloat rot[16];
float	cthru;
        if (cell->parent->empty) return; // don't draw a missing child
        if (cell->empty) return;         // or an empty shell
	chain = param->chain;
	cthru = 0.3;
	if (Data::hidden) { // distant cells become opaque
		cthru = setCthru(cell);
		if (cthru < 0.0) return;
	}
	colev = depth-level;
	col = cols[colev]; // default colour scheme
	if (moltype==0 && colours!=2) {
		if (level==depth-1) { // colour prot by SSE
			if (sort == 0) col = 3; // loop=cyan
			if (sort == 1) col = 2; // alpha=red
			if (sort == 2) col = 1; // beta=green)
		} else {	// skip over red+green+cyan
			if (level<depth) col = cols[colev+3];
		}
	}
	if (moltype==1) { // nucleic
		if (level==depth-1 && sort==0) col = 3; // loop=cyan
	}
	if (moltype==3) { // for cells, make the leading bead a different colour
		if (level==depth && cell->colour > -1) {
			col = cell->colour;
		} else {
			col = model+1;
		}
	}
	if (colours==2) { // flash bumps by cousin level: cyan < green < yellow < red < magenta < white
		col = 0; // dark blue (good to flash against)
		if (hit) { Cell *mypa, *itpa;
			if (level>0) {
				mypa = cell->parent; itpa = hit->parent;
				 if (mypa == itpa) col = 3;	// cyan
			}
			if (col==0 && level>1) {
				mypa = mypa->parent; itpa = itpa->parent;
				 if (mypa == itpa) col = 1;	// green
			}
			if (col==0 && level>2) {
				mypa = mypa->parent; itpa = itpa->parent;
				 if (mypa == itpa) col = 5;	// yellow
			}
			if (col==0 && level>3) {
				mypa = mypa->parent; itpa = itpa->parent;
				 if (mypa == itpa) col = 2;	// red
			}
			if (col==0 && level>4) {
				mypa = mypa->parent; itpa = itpa->parent;
				 if (mypa == itpa) col = 4;	// magenta
			}
			if (col==0 && level>5) {
				mypa = mypa->parent; itpa = itpa->parent;
				 if (mypa == itpa) col = 8;	// white
			}
		}
	}
	if (colours==3) { // flash Hbonds green
		if (level==depth && cell->link && cell->link[0].to != 0 && cell->link[0].type==16) col = 1;
	
	}
	if (colours==4) { // colour heat level (black--white)
		if (heat > 1.0) heat = 1.0;
		FOR(i,3) colour[i][1][0]=colour[i][1][1]=colour[i][1][2] = 1.0-heat;
		col = 1;
	}
	if (solid>0 || level==depth) {	// render solid object
		if (moltype==1 && level==depth) {
			// connect nucleic backbone (looks silly with just basepairs)
			if (parent->parent->solid > 0) return;
		} else {
			if (parent->solid > 0) return;
		}
		glMaterialfv(GL_FRONT, GL_DIFFUSE, colour[2][col]);
	} else {			  		// render cthru object
		if (parent->solid > 0) cthru *= 0.1;	// fade out contents of solid parent
		colour[0][col][3] = cthru;
		glMaterialfv(GL_FRONT, GL_DIFFUSE, colour[0][col]);
		glDepthMask(GL_FALSE);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	glPushMatrix();
		if (type==1 && chain[level]==0) { // unlinked sphere
			glTranslatef(cell->xyz.x, cell->xyz.y, cell->xyz.z);
		} else {	// draw anything with an orientation
			glTranslatef(cell->xyz.x, cell->xyz.y, cell->xyz.z);
			for (i=0; i<15; i++) rot[i] = 0.0; rot[15] = 1.0;
			rot[ 0] = cell->rot.A.x; rot[ 1] = cell->rot.A.y; rot[ 2] = cell->rot.A.z;
			rot[ 4] = cell->rot.B.x; rot[ 5] = cell->rot.B.y; rot[ 6] = cell->rot.B.z;
			rot[ 8] = cell->rot.C.x; rot[ 9] = cell->rot.C.y; rot[10] = cell->rot.C.z;
			glMultMatrixf(rot);
		}
		t = s = size;	// t = thick, s = length
		if (type==0) s = t = 1.0;	// don't scale virtual spheres
		if (type==1) rat = 0;	// spheres (just one shape type)
		if (type==2) // SSE tubes, types series = rad:len; 0 = 1:1, 1 = 1:2, ... 9 = 1:10
		{ float have, want;
			s = cell->endN | cell->endC; // length
			if (moltype==0) // PROT
			{ float	scale[3] = {loopTHIC, alphTHIC, betaTHIC};
				if (subtype==0 && level==depth-1) t *= scale[sort]; // apply SSE scale
				r = 0.5*t;		// radius
				want = s/r;		// actual len/rad ratio
				rat = (int)(0.5+want);	// rat-1 = best fit tube type
				if (rat > 10) rat = 10;
				if (rat < 1 ) rat = 1;
				have = (float)rat;	// ideal tube len/rad
				s = size * want/have;	// rescale factor to get exact tube length
				s *= scale[sort];
			}
			if (moltype==1) { // Nucleic (use given length from end points)
				r = 0.5*t;		// radius
				want = s/r;		// actual len/rad ratio
				rat = (int)(0.5+want);	// rat-1 = best fit tube type
				if (rat > 10) rat = 10;
				if (rat < 1 ) rat = 1;
				have = (float)rat;	// ideal tube len/rad
				s = size * want/have;	// rescale factor to get exact tube length
			}
			rat--;
		}
		if (type==3) // ellipsoid series = a:b; 0 = 1:10, .oblate. 10 = 1:1, .prolate. 20 = 10:1
		{ float	z = Data::Eratio[sort];		// ideal ellipsoid (half) axis length
			s = cell->endN | cell->endC;	// current axis length
			if (cell->ends < 0) t = s/z;	// the given axis length sets the size
			s = s/z;		// scale factor for exact axis length
			rat = sort;
		}
		glScalef(t,t,s); // base shape: thick = 1, length = 0.5*(rat+1)
		if (type==0) rat = 0;
		gran = log(cell->far);
		if (gran<0) gran = 0; if (gran > N-1) gran = N-1;
		glCallList(types[model][level][rat].body[gran]);
	glPopMatrix();
	if (nbonds && cell->nbonds) { // draw bonds in a chain
		if (moltype==3) col = model+1;
		if (colours==4) col = 1; // b/w heat
		if (solid>0 || level==depth || cell->live==1) {	// render solid bond
			glMaterialfv(GL_FRONT, GL_DIFFUSE, colour[2][col]);
		} else {			  		// render cthru bond
			colour[0][col][3] = cthru;
			glMaterialfv(GL_FRONT, GL_DIFFUSE, colour[0][col]);
			glDepthMask(GL_FALSE);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		}
		drawBonds(cell,cell->sis,cell->bond[0].link);
		for (int i=1; i<=nbonds; i++) { Bonds *bi = cell->bond+i;
			if (bi->to==0) continue;
			// draw hinge above 0 (sis)
			if (bi->link > 0) drawBonds(cell,bi->to,bi->link);
		}
		if (nbonds>1) { // draw crosslinks (any bond->link > 0 is a chain link)
			for (int i=2; i<=nbonds; i++) { Bonds *bi = cell->bond+i;
				if (bi->to==0) continue;
				if (bi->link < 0) drawBonds(cell,bi->to,bi->link);
			}
		}
	}
	if (Data::model[model].links[level]) drawLinks(cell);
	glDepthMask(GL_TRUE);
}

void rankLive ()
{
Cell	*world = Cell::world;
int	*p, i, n = world->kids;
int	full = 999; // 6; // additional members kept fully live
float	*far;
	p = (int*)alloca(sizeof(int)*n);
	far = (float*)alloca(sizeof(float)*n);
	for (i=0; i<n; i++) far[i] = -world->child[i]->far;
	sort(far,p,n);
	for (i=0; i<n; i++) { int live = depth - i + full;
		if (live < 1) live = 1;
		if (world->child[p[i]]->parent->type <= 0) live++; // give child locoparentis
		world->child[p[i]]->live = live;
	}
}

int getFront ()
{
float	f, x,y,z;
int	i, max,	mark[6] = {  1,  2,  3, -1, -2, -3 };
Vec	view[6];
GLfloat rot[16];
	glGetFloatv( GL_MODELVIEW_MATRIX, rot );
	view[0].x = rot[ 0]; view[0].y = rot[ 1]; view[0].z = rot[ 2];
	view[1].x = rot[ 4]; view[1].y = rot[ 5]; view[1].z = rot[ 6];
	view[2].x = rot[ 8]; view[2].y = rot[ 9]; view[2].z = rot[10];
	view[3].x = -view[0].x; view[3].y = -view[0].y; view[3].z = -view[0].z;
	view[4].x = -view[1].x; view[4].y = -view[1].y; view[4].z = -view[1].z;
	view[5].x = -view[2].x; view[5].y = -view[2].y; view[5].z = -view[2].z;
	f = -BIG; max = 0;
	for (i=0; i<6; i++) {
		if (view[i].z > f) { f=view[i].z; max=i; }
	}
	return mark[max];
}

void drawWorld ( void )
{
Cell	*world = Cell::world;
Cell	*scene = Cell::scene;
int	front, flip, use; // 0=X, 1=Y, 2=Z;
int	i,j,k, m,n;
float	d;
	if (Data::frozen < 0) { // track specified focus (slowly)
		float r = 0.01, t = 1.0-r;
		focus = focus*t + Data::focus*r;
	}
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	trX += tx; trY += ty; trZ += tz;
	rtX += rx; rtY += ry; rtZ += rz;
	tx = ty = tz = 0.0;		// one key click move
	if (Data::model[0].moltype==3) {	// if cells, fill the ring with triangles
		for (i=0; i<world->kids; i++) { Cell *cell = world->child[i];
                	if (cell->kids > 2 ) fillCell(cell);
		}
	}
	orientBody(Cell::world);	// sets the current rotation for the cell body (and link ends))
	glLoadIdentity();
	//R.set_rot('X',0.00001); // clockwise (LH) around X viewed from right of window
	//R.set_rot('Y',0.00001); // clockwise (LH) around Y viewed from top of window
	//R.set_rot('Z',0.00001); // clockwise around Z looking into window
	View *= R;		// rotate the viewing frame by R
	eye.x = trX; eye.y = trY; eye.z = trZ-dout;
	eye = View*eye;		// rotate the eye position into the viewing frame
	gluLookAt(eye.x,eye.y,eye.z, focus.x,focus.y,focus.z, View.B.x,View.B.y,View.B.z);
	rankLive();
	front = getFront();		// sets global int front (-3..+3)
	if (front > 0) { 
		flip =  1; use =  front-1;	// render from front
		scene->sort =  use+1;
	} else { 
		flip = -1; use = -front-1;	// render from back
		scene->sort = -use-1;
	}
	// current axis (+1) passed back in Cell::scene->sort for sorting in sorter()
	if (flip > 0 ) {
		for  ( i=0;  i<total;  i++ )  drawCell(scene->rank[use][i]);
	} else {
		for (i=total; i<total*2; i++) drawCell(scene->rank[use][i]);
	}
	glutSwapBuffers();
}

void spinDisplay(int value)
{
	glutPostRedisplay();
	glutTimerFunc(10,spinDisplay,0);
}

void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(30, (GLfloat) w/(GLfloat) h, 1.0, back);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity ();
	gluLookAt(0,0,dout, 0,0,0, 0,1,0);
}

void keyboard(unsigned char key, int x, int y)
{
float	r, rate = 0.0001, tate = 0.5;
	switch (key) {
		case 'X': tx += tate; break;
		case 'x': tx -= tate; break;
		case 'Y': ty += tate; break;
		case 'y': ty -= tate; break;
		case 'Z': tz += tate; break;
		case 'z': tz -= tate; break;
		// clicks inc/dec the rotation rate in R (applied to View for LookAt)
		case 'J': rx += rate; R.set_rot('X', rate); break;
		case 'j': rx -= rate; R.set_rot('X',-rate); break;
		case 'K': ry += rate; R.set_rot('Y', rate); break;
		case 'k': ry -= rate; R.set_rot('Y',-rate); break;
		case 'L': rz += rate; R.set_rot('Z', rate); break;
		case 'l': rz -= rate; R.set_rot('Z',-rate); break;
		case 'q': rx=ry=rz=0.0; hold=1;  R.ident(); break;
		case 27: printf("Viewing done\n");
			if (Data::putpdb) {
				Pt(PDB file written to final.pdb\n)
				putpdb("final.pdb");
			}
			exit(0);
	}
}

void viewer ( int argcin, char** argvin )
{
	total = Cell::total;
	depth = Data::depth;
	focus = Vec(0,0,0); // Data::focus;
	dout = 50.0; // + sqrt(total);
	glutInit(&argcin, argvin);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_ALPHA);
	glutInitWindowSize(500, 500);
	glutCreateWindow(argvin[0]);
	init();
	types = (Types***)malloc(sizeof(Types**)*N); TEST(types) // allocate for types
	makeTypes();
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutDisplayFunc(drawWorld);
	glutTimerFunc(0,spinDisplay,0);
	glutMainLoop();
}
