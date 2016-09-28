class Vec;
class Mat;
class Seg;

#include "Vec.hpp"
#include "Mat.hpp"
#include "Seg.hpp"

// vector functions
Vec   get_rand ();
Vec   get_rand ( const float );
float vdif  ( Vec, Vec );		// old style
Vec   vave ( Vec a, Vec b );
void  vave ( Vec a, Vec b, Vec *c );	// old style
float vtri  ( Vec, Vec, Vec );
float pdotp ( Vec, Vec, Vec, Vec );
float pvol  ( Vec, Vec, Vec, Vec );
float phand ( Vec, Vec, Vec, Vec );

// matrix functions
Mat frame ( Vec, Vec, Vec );
// sets an orthogonal coordinate frame on b with y as the bisector of abc and x close to ac

// segment functions

float dist_to_cut ( const Vec&, const Vec&, const Vec& );
// returns the distance to the image of c on the extended line along a->b (as a fraction of |a-b|)

int  parallel ( Vec, Vec );
// return 1 if vectors a and b are parallel and -1 if antiparallel

int  parallel ( Seg, Seg );
// return -1/1 if line segments a and b are anti/parallel

bool orthgnal ( Vec, Vec );
// return 1 if vectors a and b are orthogonal (within +/-NOISE)

bool orthgnal ( Seg, Seg );
// return 1 if line segments a and b are orthogonal

Vec getz ( Vec, Vec, Vec, Vec );
// returns the normal vector between the lines at the middle of their overlap segment 

Vec get_comp ( Vec, Vec, Vec, Vec );
// return the components of the basis set (x,y,z) to locate the contact normal.

bool line2line ( Seg, Seg, float ); // as in old geom.c
// TRUE if line segments a and b overlap and are closer than s 

bool seg_in_seg ( Seg, Seg );
// TRUE if line segments a and b overlap (ie an end has a normal in the other segment)

float line_to_line ( Seg, Seg );
// returns the closest approach distance between two extended line segments

Seg line_on_line ( Seg, Seg );
// returns the endpoints of the contact normal between extended a and b

float seg_to_seg ( Seg, Seg ); // was dist2line() in old geom.c
// returns the closest approach distance between two line segments

Seg seg_on_seg ( Seg, Seg );
// returns the endpoints of the shortest line between the two segments

// angle functions

float angle ( Vec, Vec, Vec );
// returns the (unsigned 0...PI) angle at b (in rad)

float angle_deg ( Vec, Vec, Vec );
// returns the (unsigned 0...180) angle at b (in deg)

float angle1pi ( float, float );
// returns the angle (-PI...0...PI) from sin and cos components

float angle2pi ( float, float );
// returns the angle (0...2PI) from sin and cos components

float angdif ( float, float );
// absolute difference of two angles (in +/- PI range)

float torsion ( Vec, Vec, Vec, Vec );
// returns the torsion angle down b-c given a,b,c,d

// separate function
void separate ( Vec&, Vec&, float, float ); // pass by reference
void separate ( Vec*, Vec*, float, float ); // pass by pointer

// ellips[e/oid] functions

float toEllipse ( float, float, float, float );
// distance of a point to the perimiter of an ellipse (in 2D)

float vec_to_egg ( Vec, Seg, float );
float vec_to_egg ( Vec, Vec, Vec, float );
// distance of a point to the perimiter of an ellipse (in 3D)
