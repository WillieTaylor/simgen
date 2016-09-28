#include "util.hpp"
#include "geom.hpp"

// vector functions
void Vec::set_disp () {                      // add a random displacement between +/-1 to each component
        x += 2.0*(randf()-0.5); y += 2.0*(randf()-0.5); z += 2.0*(randf()-0.5);
}
void Vec::set_disp ( const float s ) {       // add a random displacement between +/-<s> (was vradd())
        x += s*2.0*(randf()-0.5); y += s*2.0*(randf()-0.5); z += s*2.0*(randf()-0.5);
}

void Vec::set_rand () {			// set a random vector between +/-1
        x = 2.0*(randf()-0.5); y = 2.0*(randf()-0.5); z = 2.0*(randf()-0.5);
}
void Vec::set_rand ( const float s ) {	// set a random vector between +/-<s>
        x = s*2.0*(randf()-0.5); y = s*2.0*(randf()-0.5); z = s*2.0*(randf()-0.5);
}

void Vec::set_frac ( const Mat &M ) {
// returns the fractional components of <this> in the basis set <M> (M need not be orthogonal)
Vec	p = (*this);
float   d;
Mat	W;
        d = M.det();
        if (fabs(d) < NOISE) { this->zero(); return; } // zero = fail
        W = M.get_inv(d);
        (*this) = W*p;
}

Vec Vec::get_frac ( const Mat &M ) {
// returns the fractional components of <p> in the basis set <M> (M need not be orthogonal)
Vec	q = (*this);
        q.set_frac(M);
	return q;
}

Vec Vec::get_frac ( const Vec &p, const Mat &M ) {
// returns the fractional components of <p> in the basis set <M> (M need not be orthogonal)
Vec	q = p;
        q.set_frac(M);
	return q;
}

void Vec::set_frac ( const Vec &p, const Mat &M ) {
// returns the fractional components of <this> in the basis set <M> (M need not be orthogonal)
	(*this) = get_frac(p,M);
}

Vec get_rand () {		// retrun a random displacement between +/-1
        return Vec ( 2.0*(randf()-0.5), 2.0*(randf()-0.5), 2.0*(randf()-0.5) );
}
Vec get_rand ( const float s ) { // retrun a random displacement between +/-<s>
        return Vec ( s*2.0*(randf()-0.5), s*2.0*(randf()-0.5), s*2.0*(randf()-0.5) );
}

float vdif ( Vec a, Vec b ) { return a|b; }

Vec   vave ( Vec a, Vec b ) { return a&b; }
void  vave ( Vec a, Vec b, Vec *c ) { *c =  a&b; } // old style

float vtri ( Vec a, Vec b, Vec c ) { return (a^b)*c; }

float pdotp ( Vec a, Vec b, Vec c, Vec d ) { return (b-a)*(d-c); }	 

float pvol ( Vec a, Vec b, Vec c, Vec d )
{	Vec	x, y, z;
	x = b-a; y = d-c; z = c-b;
	return vtri(x,y,z);
}	 

float phand ( Vec a, Vec b, Vec c, Vec d )
{	Vec	x, y, z;
	x.setVec(b-a); y.setVec(d-c); z.setVec(c-b);
	return vtri(x,y,z);
}	 

// matrix operators
Mat Mat::operator * ( const Mat &Q ) {		// multiply two matrices
	Mat P = *this;
	A.x = P.A.x*Q.A.x + P.B.x*Q.A.y + P.C.x*Q.A.z;
       	B.x = P.A.x*Q.B.x + P.B.x*Q.B.y + P.C.x*Q.B.z;
       	C.x = P.A.x*Q.C.x + P.B.x*Q.C.y + P.C.x*Q.C.z;

	A.y = P.A.y*Q.A.x + P.B.y*Q.A.y + P.C.y*Q.A.z;
       	B.y = P.A.y*Q.B.x + P.B.y*Q.B.y + P.C.y*Q.B.z;
       	C.y = P.A.y*Q.C.x + P.B.y*Q.C.y + P.C.y*Q.C.z;

	A.z = P.A.z*Q.A.x + P.B.z*Q.A.y + P.C.z*Q.A.z;
       	B.z = P.A.z*Q.B.x + P.B.z*Q.B.y + P.C.z*Q.B.z;
       	C.z = P.A.z*Q.C.x + P.B.z*Q.C.y + P.C.z*Q.C.z;
	return Mat(A,B,C);
}
Mat Mat::operator *= ( const Mat &Q ) {		// multiply and set matrices
	Mat P = (*this) * Q;
	A = P.A; B = P.B; C = P.C;
}
Vec Mat::operator * ( const Vec &d ) {		// multiply a Mat and a Vec (was MmulV)
	Vec a = A, b = B, c = C;
	Vec e;
        e.x = d.x*a.x + d.y*b.x + d.z*c.x;
        e.y = d.x*a.y + d.y*b.y + d.z*c.y;
        e.z = d.x*a.z + d.y*b.z + d.z*c.z;
	return e;
}
// extra Vec operator
Vec Vec::operator * ( const Mat &M ) {		// multiply a Vec and a Mat (was VmulM)
	Vec a = M.A, b = M.B, c = M.C;
	Vec e, d = *this;
	e.x = d.x*a.x + d.y*a.y + d.z*a.z;
	e.y = d.x*b.x + d.y*b.y + d.z*b.z;
	e.z = d.x*c.x + d.y*c.y + d.z*c.z;
	return e;
}
Vec Vec::operator *= ( const Mat &M ) {		// multiply a Vec and a Mat and set
	Vec v = (*this) * M;
	x = v.x; y = v.y; z = v.z;
}

// matrix inverse functions
void Mat::set_inv ( const float det ) {		// set the current matrix to its inverse (given det)
	Vec a = A, b = B, c = C;
	if (fabs(det)<NOISE) {
		cout << "Singular matrix in set_inv()\n";
		(*this).print();
		return;
	}
	A.x =  (b.y*c.z - b.z*c.y)/det;
	A.y = -(a.y*c.z - a.z*c.y)/det;
	A.z =  (a.y*b.z - a.z*b.y)/det;

	B.x = -(b.x*c.z - b.z*c.x)/det;
	B.y =  (a.x*c.z - a.z*c.x)/det;
	B.z = -(a.x*b.z - a.z*b.x)/det;

	C.x =  (b.x*c.y - b.y*c.x)/det;
	C.y = -(a.x*c.y - a.y*c.x)/det;
	C.z =  (a.x*b.y - a.y*b.x)/det;
}
void Mat::set_inv () {				// set the current matrix to its inverse (find det)
	float det = (*this).det();
	(*this).set_inv(det);
}

Mat Mat::get_inv (const float det ) const {	// return the inverse of the current matrix (given det)
	Vec a = A, b = B, c = C; Mat W;
	if (fabs(det)<NOISE) {
		cout << "Singular matrix in get_inv()\n";
		(*this).print();
		return *this;
	}
	W.A.x =  (b.y*c.z - b.z*c.y)/det;
	W.A.y = -(a.y*c.z - a.z*c.y)/det;
	W.A.z =  (a.y*b.z - a.z*b.y)/det;

	W.B.x = -(b.x*c.z - b.z*c.x)/det;
	W.B.y =  (a.x*c.z - a.z*c.x)/det;
	W.B.z = -(a.x*b.z - a.z*b.x)/det;

	W.C.x =  (b.x*c.y - b.y*c.x)/det;
	W.C.y = -(a.x*c.y - a.y*c.x)/det;
	W.C.z =  (a.x*b.y - a.y*b.x)/det;
	return W;
}
Mat Mat::get_inv () const {			// return the inverse of the current matrix (find det)
	float det = (*this).det();
	return (*this).get_inv(det);
}
Mat Mat::inv ( const float det ) const { return (*this).get_inv(det); }
Mat Mat::inv () const { return (*this).get_inv(); }

void Mat::set_rot ( const char axis, const float r )
{				// rotates the current matrix by <r> rad.s about <axis> = ['X'|'Y'|'Z']
float	cr = cos(r), sr = sin(r);
Vec	a,b,c;
Mat	R;
	if (axis=='X') {   a.x = 1.0;
		b.y =  cr; b.z = -sr;
		c.y =  sr; c.z =  cr;
	}
	if (axis=='Y') {   b.y = 1.0;
		a.x =  cr; a.z = -sr;
		c.x =  sr; c.z =  cr;
	}
	if (axis=='Z') {   c.z = 1.0;
		a.x =  cr; a.y = -sr;
		b.x =  sr; b.y =  cr;
	}

	R.A.x = a.x*A.x + b.x*A.y + c.x*A.z;
        R.B.x = a.x*B.x + b.x*B.y + c.x*B.z;
        R.C.x = a.x*C.x + b.x*C.y + c.x*C.z;

	R.A.y = a.y*A.x + b.y*A.y + c.y*A.z;
        R.B.y = a.y*B.x + b.y*B.y + c.y*B.z;
        R.C.y = a.y*C.x + b.y*C.y + c.y*C.z;

	R.A.z = a.z*A.x + b.z*A.y + c.z*A.z;
        R.B.z = a.z*B.x + b.z*B.y + c.z*B.z;
        R.C.z = a.z*C.x + b.z*C.y + c.z*C.z;

	A = R.A; B = R.B; C = R.C;
}

Mat frame ( Vec a, Vec b, Vec c ) {
// sets an orthogonal coordinate frame on b with y as the bisector of abc and x close to ac
Vec	x = c-a, y = b-(c&a), z = x^y;
float	d = Mat(a,b,c).det();
	if (fabs(d)<NOISE) {
		cout << "Bad vectors for frame():  "; Pv(a) Pv(b) Pv(c) NL;
	}
	x = y^z;
	return Mat( x.norm(), y.norm(), z.norm() );
}

// segment functions

float dist_to_cut ( const Vec &a, const Vec &b, const Vec &c ) {
// returns the distance to the image of c on the extended line along a->b (as a fraction of |a-b|)
// p*q = |p||q|cos(a), D = |q|cos(a), = p*q/|p|, d = D/|p| = p*q/|p|^2 
Vec	p = b-a, q = c-a;
float	d = (p*q)/p.sqr();
	return d;
}

bool Seg::vec_in_seg ( const Vec &c ) const {
// TRUE if c lies over the line segment
float	d = dist_to_cut(A,B,c);
	if ( d<0.0 || d>1.0 ) return 0; else return 1;
}
bool Vec::vec_in_seg ( const Seg &s ) const {
// TRUE if c lies over the line segment
float	d = dist_to_cut(s.A,s.B,*this);
	if ( d<0.0 || d>1.0 ) return 0; else return 1;
}
bool Vec::vec_in_seg ( const Vec &a, const Vec &b ) const {
	return this->vec_in_seg(Seg(a,b));
}

float Seg::vec_to_line ( const Vec &c ) const {
// returns the distance from point c to the extended line segment (A-B)
Vec	p = B-A, q = c-A;
float	d = (p*q)/p.len(), // NB not sqr (as d not fractional)
	h = sqrt(q.sqr()-d*d);
	return h;
}
float Vec::vec_to_line ( const Seg &s ) const {
// returns the distance from <this> point to the extended line segment (A-B)
Vec	p = s.B-s.A, q = (*this)-s.A;
float	d = (p*q)/p.len(), // NB not sqr()
	h = q.sqr()-d*d;
	if (h < NOISE) return 0.0;
	return sqrt(h);
}
float Vec::vec_to_line ( const Vec &a ) const {
// returns the distance from <this> point to the extended line segment (0-a)
	return this->vec_to_line(Seg(Vec(0,0,0),a));
}
float Vec::vec_to_line ( const Vec &a, const Vec &b ) const {
// returns the distance from <this> point to the extended line segment (a-b)
	return this->vec_to_line(Seg(a,b));
}

Vec Seg::vec_on_line ( const Vec &c ) const { // was dot2line()
// returns the image of c on extended A-B line segment
Vec	p = B-A, q = c-A;
float	d = (p*q)/p.sqr();
Vec	g = A + p*d;
	return g;
}
Vec Vec::vec_on_line ( const Seg &s ) const { // was dot2line()
// returns the image of c on extended A-B line segment
Vec	p = s.B-s.A, q = (*this)-s.A;
float	d = (p*q)/p.sqr();
Vec	g = s.A + p*d;
	return g;
}
Vec Vec::vec_on_line ( const Vec &a, const Vec &b ) const { // as above
	return this->vec_on_line(Seg(a,b));
}

Vec Seg::hit_trif ( const Vec &a, const Vec &b, const Vec &c ) const 
{ // returns the fractional basis vector components for the intersection with the plane a,b,c (was part of line2tri())
Vec	x = b-a, y =c-a, z = A-B;
Mat	M = Mat(x,y,z);
float	det = M.det();
Vec	d,e;
	if (fabs(det) < NOISE) return Vec(); // zero = fail
	M.set_inv(det);
	d = A-a; e = M*d;
	return e;
}

Vec Seg::hit_tri ( const Vec &a, const Vec &b, const Vec &c ) const 
{ // returns the intersection point of the extended line segment and the plane a,b,c
	Vec d = Vec( b|a, c|a, A|B);
	Vec e = hit_trif(a,b,c);
	return Vec( a.x+d.x*e.x, a.y+d.y*e.y, a.z+d.z*e.z );
}

Vec Seg::hit_tri ( const Vec &a, const Vec &b, const Vec &c, const Vec &d ) const 
{ // returns the intersection point of the extended line segment and the plane a,b,c (given d)
	Vec e = Vec( b|a, c|a, A|B);
	return Vec( a.x+d.x*e.x, a.y+d.y*e.y, a.z+d.z*e.z );
}

bool Seg::cut_tri ( const Vec &a, const Vec &b, const Vec &c ) const 
{	// TRUE if segment intersects the triangle a,b,c (was line2tri())
	Vec e = hit_trif(a,b,c);
	if (e.x < 0.0) return 0;
	if (e.y < 0.0) return 0;
	if (e.z < 0.0) return 0;
	if (e.z > 1.0) return 0;
	if (e.x + e.y > 1.0) return 0;
	return 1;
}

bool Seg::cut_tri ( const Vec &e ) const 
{	// TRUE if segment intersects the triangle a,b,c (was line2tri())
	if (e.x < 0.0) return 0;
	if (e.y < 0.0) return 0;
	if (e.z < 0.0) return 0;
	if (e.z > 1.0) return 0;
	if (e.x + e.y > 1.0) return 0;
	return 1;
}

int parallel ( Vec x, Vec y ) {
// return 1 if vectors a and b are parallel and -1 if antiparallel
float	cos = x.norm()*y.norm();
	if (cos >  1.0-NOISE) return  1;
	if (cos < -1.0+NOISE) return -1;
	return 0;
}

int parallel ( Seg a, Seg b ) {
// return -1/1 if line segments a and b are anti/parallel
Vec	x = a.vec(), y = b.vec();
	return parallel(x,y);
}

bool orthgnal ( Vec x, Vec y ) {
// return 1 if vectors a and b are orthogonal
float	cos = x.norm()*y.norm();
	if (cos > -NOISE && cos < NOISE) return  1;
	return 0;
}

bool orthgnal ( Seg a, Seg b ) {
// return 1 if line segments a and b are orthogonal
Vec	x = a.vec(), y = b.vec();
	return orthgnal(x,y);
}

Vec getz ( Vec a1, Vec b1, Vec a2, Vec b2 ) {
// for parallel lines a1-->b1, a2-->b2, 
// returns the normal vector between the lines at the middle of their overlap segment 
// returns 0,0,0 if no overlap
float	ca, cb, ta, tb, d1, d2, x1, x2, h1, h2;
Vec	x, y, z, z1, z2;
	x.setVec(b1-a1); y.setVec(b2-a1); z.setVec(b1-a2);
	ca = x*y, cb = x*z;
	if ( ca < 0.0 || cb < 0.0 ) return Vec();	//  an obtuse angle = no overlap
	// tan1 = h/x, tan2 = h/(d-x), tan1/x = tan2/(d-x), tan2/tan1 = (d-x)/x = d/x - 1
	// d/x = tan2/tan1 + 1, x = d/(tan2/tan1 + 1)
	// for seg 1
	ta = tan(acos(ca)); tb = tan(acos(cb)); d1 = a1|b1;
	x1 = d1/(ta/tb + 1); h1 = ta/x1;
	z1 = a1 + x*x1;
	// for seg 2
	x.setVec(b2-a2); y.setVec(b1-a2); z.setVec(b2-a1);
	ca = x*y, cb = x*z;
	ta = tan(acos(ca)); tb = tan(acos(cb)); d2 = a2|b2;
	x2 = d2/(ta/tb + 1); h2 = ta/x2;
	z2 = a2 + x*x2;
	return z2-z1;
}

Vec	get_comp ( Vec x, Vec y, Vec z, Vec d ) {
// return the components of the basis set (x,y,z) to locate the contact normal.
// Find the components that connect the endpoints (d) through the mutual perpendicular (z)
// z = x^y. The path to connect the ends d = p*x+q*z+r*y is true in each dimension so if 
// e = {p,q,r} then d = e*M  and e can be solved for simultaneously as e = W*d (W = M^-1)
// By making z a unit vector, its component (e.z) will be the distance between the lines.
Mat	M = Mat(x,y,z.norm()),
	W = M.inv();
	return W*d;
}

bool line2line ( Seg a, Seg b, float s ) { // as in old geom.c
// TRUE if line segments a and b overlap and are closer than s 
Vec	x = a.vec(), y = b.vec(), z, d = b.B-a.A, e;
int	p = parallel(x,y);
	if (p > 0) z = getz(a.A,a.B,b.A,b.B);
	if (p < 0) z = getz(a.A,a.B,b.B,b.A);
	if (p == 0) z = x^y;
	e = get_comp(x,y,z,d);
        if (e.x < 0.0) return 0;
        if (e.y < 0.0) return 0;
        if (e.x > 1.0) return 0;
        if (e.y > 1.0) return 0;
        if (e.z >  s ) return 0;
        if (e.z < -s ) return 0;
        return 1;
}

bool seg_in_seg ( Seg a, Seg b ) {
// TRUE if line segments a and b overlap (ie an end has a normal in the other segment)
	if (a.vec_in_seg(b.A)) return 1;
	if (a.vec_in_seg(b.B)) return 1;
	if (b.vec_in_seg(a.A)) return 1;
	if (b.vec_in_seg(a.B)) return 1;
	return 0;
}

float line_to_line ( Seg a, Seg b ) {
// returns the closest approach distance between two extended line segments
Vec	x = a.vec(), y = b.vec(), z, d = b.B-a.A, e;
	if (parallel(x,y)) { // return the distance of any point to the other line
		return (a.A).vec_to_line(b);
	}
	z = x^y;
	e = get_comp(x,y,z,d);
	return fabs(e.z);
}

Seg line_on_line ( Seg a, Seg b ) {
// returns the endpoints of the contact normal between extended a and b
// for parallel lines, the normal is in the middle of the overlap region
// or if no overlap, passes through the centroid of the 4 endpoints
Vec	x = a.vec(), y = b.vec(), z, d = b.B-a.A, e;
int	p = parallel(x,y);
	if (p > 0) z = getz(a.A,a.B,b.A,b.B);
	if (p < 0) z = getz(a.A,a.B,b.B,b.A);
	if (p && z.iszero()) { Vec mid = a.mid() & b.mid(); // on overlap
		return Seg( mid.vec_to_line(a), mid.vec_to_line(b) );
	}
	if (p == 0) z = x^y;
	e = get_comp(x,y,z,d);
	return Seg( a.A+x*e.x, b.B-y*e.y );
}

float seg_to_seg ( Seg a, Seg b ) { // was dist2line() in old geom.c
// returns the closest approach distance between two line segments
float	d, dmin = BIG, dend;
Vec	x = a.vec(), y = b.vec(), z, f = b.B-a.A, e;
int	p = parallel(x,y);
	if (p > 0) z = getz(a.A,a.B,b.A,b.B);
	if (p < 0) z = getz(a.A,a.B,b.B,b.A);
	if (p && z.iszero()) { // on overlap, return closest pair separation
		return fmin(fmin(a.A|b.A,a.A|b.B),fmin(a.B|b.A,a.B|b.B));
	}
	if (p==0) z = x^y;
	e = get_comp(x,y,z,f);
        if ((e.x>0.0 && e.x<1.0) && (e.y>0.0 && e.y<1.0)) { // normal lies in both segments
		// since z has unit length (in get_comp()) its component is the real length
		return fabs(e.z);
	}
	// find the closest end ([ab].[AB] over a line [ab]
	if (a.vec_in_seg(b.A)) dmin = fmin(dmin,a.vec_to_line(b.A));
	if (a.vec_in_seg(b.B)) dmin = fmin(dmin,a.vec_to_line(b.B));
	if (b.vec_in_seg(a.A)) dmin = fmin(dmin,b.vec_to_line(a.A));
	if (b.vec_in_seg(a.B)) dmin = fmin(dmin,b.vec_to_line(a.B));
	if (dmin < BIG-1.0) return dmin;
	// otherwise find the closest ends
	dend = fmin(fmin(a.A|b.A,a.A|b.B),fmin(a.B|b.A,a.B|b.B));
	if (dend < dmin) dmin = dend;
	return dmin;
}

Seg seg_on_seg ( Seg a, Seg b ) {
// returns the endpoints of the shortest line between the two segments
Seg	best;
float	d, dmin = BIG;
Vec	x = a.vec(), y = b.vec(), z, f = b.B-a.A, e;
int	p = parallel(x,y);
	if (p > 0) z = getz(a.A,a.B,b.A,b.B);
	if (p < 0) z = getz(a.A,a.B,b.B,b.A);
	if (p && z.iszero()) { // on overlap, return closest pair separation
		d = a.A|b.A; if (d<dmin) { dmin=d; best.A=a.A; best.B=b.A; }
		d = a.A|b.B; if (d<dmin) { dmin=d; best.A=a.A; best.B=b.B; }
		d = a.B|b.A; if (d<dmin) { dmin=d; best.A=a.B; best.B=b.A; }
		d = a.B|b.B; if (d<dmin) { dmin=d; best.A=a.B; best.B=b.B; }
		return best;
	}
	if (p==0) z = x^y;
	e = get_comp(x,y,z,f);
        if ((e.x>0.0 && e.x<1.0) && (e.y>0.0 && e.y<1.0)) { // normal lies in both segments
		best.A = a.A+x*e.x; best.B = b.B-y*e.y;
		return best;
	}
	// find the closest end over a line
	if (a.vec_in_seg(b.A)) { d = a.vec_to_line(b.A); if (d<dmin) { dmin=d; best.A=a.vec_on_line(b.A); best.B=b.A; } }
	if (a.vec_in_seg(b.B)) { d = a.vec_to_line(b.B); if (d<dmin) { dmin=d; best.A=a.vec_on_line(b.B); best.B=b.B; } }
	if (b.vec_in_seg(a.A)) { d = b.vec_to_line(a.A); if (d<dmin) { dmin=d; best.B=b.vec_on_line(a.A); best.A=a.A; } }
	if (b.vec_in_seg(a.B)) { d = b.vec_to_line(a.B); if (d<dmin) { dmin=d; best.B=b.vec_on_line(a.B); best.A=a.B; } }
	if (dmin < BIG-1.0) return best;
	// otherwise find the closest ends
	d = a.A|b.A; if (d<dmin) { dmin=d; best.A=a.A; best.B=b.A; }
	d = a.A|b.B; if (d<dmin) { dmin=d; best.A=a.A; best.B=b.B; }
	d = a.B|b.A; if (d<dmin) { dmin=d; best.A=a.B; best.B=b.A; }
	d = a.B|b.B; if (d<dmin) { dmin=d; best.A=a.B; best.B=b.B; }
	return best;
}

float Seg::operator  | ( const Seg &s ) const {
// if the two segments overlap return their closest approach, otherwise -1
	if (seg_in_seg(*this,s)) return seg_to_seg(*this,s); else return -1.0;
}
float Seg::operator || ( const Seg &s ) const {	// squared distance
	if (seg_in_seg(*this,s))
	{ float d = seg_to_seg(*this,s); return d*d;
	} else { return -1.0; }
}

Vec Vec::get_rot ( const Seg &s, const float t ) const {
// return current point (this) rotated about a-b by t (+t = clockwise viewed along a-->b)
Vec	c = *this, d = c-s.A;
Vec	x = s.B-s.A, y = d, z;
Mat	M,R = Mat();
	R.B.y = cos(t); R.B.z = -sin(t);
	R.C.y = sin(t); R.C.z =  cos(t);
	if (x.iszero()) return c; // b lies on a
	if (y.iszero()) return c; // c lies on a
	z = x^y;
	if (z.iszero()) return c; // c lies on a--b
	y = x^z;
	M = Mat(x.norm(),y.norm(),z.norm());
	d = d*M;		// a-b to X
	d = d*R;		// rotate around X by t
	d = M*d;		// X to a-b (post-mul = transpose = inverse)
	d += s.A;
	return d;
}
Vec Vec::get_rot ( const Vec &a, const Vec &b, const float t ) const {
	return this->get_rot(Seg(a,b),t);
}

void Vec::set_rot ( const Seg &s, const float t ) {
// set current point (this) rotated about a-b by t (+t = clockwise viewed along a-->b)
Vec	c = this->get_rot(s,t);
	x = c.x; y=c.y; z=c.z;
}
void Vec::set_rot ( const Vec &a, const Vec &b, const float t ) { // dito using two vectors
Vec	c = this->get_rot(Seg(a,b),t);
	x = c.x; y=c.y; z=c.z;
}
void Vec::set_rot ( const Vec &a, const float t ) { // dito using the origin and a vector
Vec	c = this->get_rot(Seg(Vec(),a),t);
	x = c.x; y=c.y; z=c.z;
}
void Vec::set_rot ( const float t ) { // dito using the origin and a random vector
Vec	c = this->get_rot(Seg(Vec(),Vec(1.0)),t);
	x = c.x; y=c.y; z=c.z;
}

// angle functions

float angle ( Vec a, Vec b, Vec c ) {
// returns the (unsigned 0...PI) angle at b (in rad)
float	d;
Vec	x = a-b, y = c-b;
	if (x.iszero()) return 0.0; // b lies on a
	if (y.iszero()) return 0.0; // b lies on c
	d = x.norm() * y.norm();
	if (d < NOISE-1.0) return PI;
	return acos(d);
}

float angle_deg ( Vec a, Vec b, Vec c ) {
// returns the (unsigned 0...180) angle at b (in deg)
	return 180.0*angle(a,b,c)/PI;
}

float angle1pi ( float s, float c ) {
// returns the angle (-PI...0...PI) from sin and cos components
	if (s>=0.0 && c>=0.0) {
		if (s<0.5) return asin(s);
		      else return acos(c);
	}
	if (s>=0.0 && c<=0.0) {
		if (s<0.5) return PI - asin(s);
		      else return PI - acos(-c);
	}
	if (s<=0.0 && c<=0.0) {
		if (s>-.5) return asin(-s) - PI;
		      else return acos(-c) - PI;
	}
	if (s<=0.0 && c>=0.0) {
		if (s>-.5) return -asin(-s);
		      else return -acos(c);
	}
	printf("angle1pi(s,c) out of range: s = %f, c = %f\n", s,c);
	return 999.999;
	//exit(1);
}

float	angle2pi ( float s, float c ) {
// returns the angle (0...2PI) from sin and cos components
	if (s>=0.0 && c>=0.0) {
		if (s<0.5) return asin(s);
		      else return acos(c);
	}
	if (s>=0.0 && c<=0.0) {
		if (s<0.5) return PI - asin(s);
		      else return PI - acos(-c);
	}
	if (s<=0.0 && c<=0.0) {
		if (s>-.5) return PI + asin(-s);
		      else return PI + acos(-c);
	}
	if (s<=0.0 && c>=0.0) {
		if (s>-.5) return twoPI - asin(-s);
		      else return twoPI - acos(c);
	}
	printf("angle2pi(s,c) out of range: s = %f, c = %f\n", s,c);
	return 999.999;
	//exit(1);
}

float angdif ( float a, float b ) {
// absolute difference of two angles (in +/- PI range)
        if (a<0.0 && b>0.0)
        { float d = b-a;
                if (d<PI) return d;
                     else return twoPI-d;
        }
        if (a>0.0 && b<0.0)
        { float d = a-b;
                if (d<PI) return d;
                     else return twoPI-d;
        }
        if (a>b) return a-b; else return b-a;
}

float torsion ( Vec a, Vec b, Vec c, Vec d ) {
// returns the torsion angle down b-c 
Vec	x = a-b, y = b-c, z = c-d, p, q, r;
float	vol, cos, sin, tor;
	if ((a|d) < NOISE) return PI;
	tor = 2.0*(randf()-0.5)*PI; // default random
	if (x.iszero()) return tor; // b lies on a
	if (y.iszero()) return tor; // b lies on c
	if (z.iszero()) return tor; // d lies on c
	p = x^y; p.setVec();
	if (p.iszero()) return tor; // x,y colinear
	q = y^z; q.setVec();
	if (q.iszero()) return tor; // y,z colinear
	cos = p*q; r = p^q;
	vol = vtri(p,q,y);
	sin = r.len();
	if (vol < 0.0) sin = -sin;
	tor = angle1pi(sin,cos);
	return tor;
}

// separate functions

void Seg::separate ( float dist, float kick )
{
Vec     oldm = A, oldn = B, disp, mid;
float   gap, dif;
        gap = oldm|oldn;
	if (gap<NOISE) { oldm = Vec(NOISE); oldn = Vec(NOISE); }
	if (kick > 1.0) kick = 1.0; // better to use bigger dist, not kick
	if (kick < 0.0) {	// -ve kick repels only
		if (gap > dist) return;
		kick = -kick;
	}
        dif = 0.5*(gap-dist);
	if (fabs(dif) < NOISE) return;
	mid = oldm & oldn;
	disp.setVec(oldm-mid);
	disp *= kick*dif;
        A = oldm-disp; B = oldn+disp;
}
void Seg::separate ( float dist ) {
	separate(dist,1.0);
}

void separate ( Vec &a, Vec &b, float dist, float kick ) {
Seg	temp = Seg(a,b);
	temp.separate(dist,kick);
	a = temp.A; b = temp.B;
}
void separate ( Vec *a, Vec *b, float dist, float kick ) {
Seg	temp = Seg(*a,*b);
	temp.separate(dist,kick);
	*a = temp.A; *b = temp.B;
}

float toEllipse ( float x, float y, float a, float b ) {
// returns the distance of 2D point <p> to an ellipse with axes A=ab.x, B=ab.y
// see: http://www.iquilezles.org/www/articles/ellipsedist/ellipsedist.htm
Vec p = Vec(x,y,0), ab = Vec(a,b,0);
    if (fabs(ab.x-ab.y) < NOISE) return p.len()-ab.x;	// sphere
    if (p.y < NOISE) return p.x-ab.x;		// on A axis
    if (p.x < NOISE) return p.y-ab.y;		// on B axis
    double ax = p.x/ab.x, by = p.y/ab.y;
    if (ax*ax + by*by < 1.0) return -1.0;	// inside
    if (p.x < 0.0) p.x = -p.x;
    if (p.y < 0.0) p.y = -p.y;
    //  input OK
    double l = ab.y*ab.y - ab.x*ab.x;
    double m = ab.x*p.x/l; double m2 = m*m;
    double n = ab.y*p.y/l; double n2 = n*n;
    double c = (m2 + n2 - 1.0)/3.0; double c3 = c*c*c;
    double q = c3 + m2*n2*2.0;
    double d = c3 + m2*n2;
    double g = m + m*n2;
    double co;
    if( d<0.0 ) {
        double p = acos(q/c3)/3.0;
        double s = cos(p);
        double t = sin(p)*sqrt(3.0);
        double rx = sqrt( -c*(s + t + 2.0) + m2 );
        double ry = sqrt( -c*(s - t + 2.0) + m2 );
        co = ( ry + sign(l)*rx + fabs(g)/(rx*ry) - m)/2.0;
    } else {
        double h = 2.0*m*n*sqrt( d );
        double s = sign(q+h)*pow( fabs(q+h), 1.0/3.0 );
        double u = sign(q-h)*pow( fabs(q-h), 1.0/3.0 );
        double rx = -s - u - c*4.0 + 2.0*m2;
        double ry = (s - u)*sqrt(3.0);
        double rm = sqrt( rx*rx + ry*ry );
        double p = ry/sqrt(rm-rx);
        co = (p + 2.0*g/rm - m)/2.0;
    }
    double si = sqrt( 1.0 - co*co );
    Vec closestPoint = Vec( ab.x*co, ab.y*si, 0 );
    return (closestPoint - p).len() * sign(p.y - closestPoint.y);
}

float vec_to_egg ( Vec p, Seg axis, float size ) {
// returns the distance of a point <p> to an ellipsoid (A=B=size, C=|axis|)
Vec	q = p.vec_on_line(axis);	// q = image of p on the axis
float	x = q|(axis.A & axis.B),	// q--centre distance
	y = p|q,			// p--axis distance
	a = 0.5*axis.len(),		// semi-axis length
	b = 0.5*size;			// other axis length
	return toEllipse(x,y,a,b);
}
float vec_to_egg ( Vec p, Vec A, Vec B, float size ) {
Seg	axis = Seg(A,B);
	return vec_to_egg(p,axis,size);
}

/* not converted

float	vdad ( Vec a, Vec b )
{	Vec	c;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	c.z = a.z + b.z;
	return sqrt(c.x*c.x + c.y*c.y + c.z*c.z);
}	 

float	vddad ( Vec a, Vec b )
{	Vec	c;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	c.z = a.z + b.z;
	return (c.x*c.x + c.y*c.y + c.z*c.z);
}	 

void	vat0 ( Vec a, Vec b, Vec *c, float d )
// c = vector of length d in direction a-->ba (at the origin)
{
Vec	e;
float	f;
	e.x = b.x - a.x;
	e.y = b.y - a.y;
	e.z = b.z - a.z;
	f = sqrt(e.x*e.x + e.y*e.y + e.z*e.z);
	if (f < NOISE) {
		c->x = c->y = c->z = 0.0;
		return;
	}
	c->x = d*e.x/f;
	c->y = d*e.y/f;
	c->z = d*e.z/f;
}

void	vatA ( Vec a, Vec b, Vec *c, float d )
// c = vector of length d in direction a-->b at a
{
Vec	e;
float	f;
	e.x = b.x - a.x;
	e.y = b.y - a.y;
	e.z = b.z - a.z;
	f = sqrt(e.x*e.x + e.y*e.y + e.z*e.z);
	if (f < NOISE) {
		c->x=a.x; c->y=a.y; c->z=a.z;
		return;
	}
	c->x = a.x + d*e.x/f;
	c->y = a.y + d*e.y/f;
	c->z = a.z + d*e.z/f;
}

int norm2line ( Vec a, Vec b, Vec c, Vec d, Vec *g, Vec *h ) {
// g,h = closest points on line segments a-b and c-d
// returns: 0 = both in line, 1..4 = end a..d, -1..-4 = a+c, a+d, b+c, b+d
Vec x, y, z;
Mat M[1], W[1];
float det, dmin;
Vec e,f;
        vsub(b,a,&x);
        vsub(d,c,&y);
        vprod(x,y,&z);
        vnorm(&z);
        VtoM(x,y,z,M);
        det = Mdet(M);
        Minv(M,W,det);
        vsub(d,a,&f);
        MmulV(W,f,&e);
	vmul(&x,e.x);
	vadd(a,x,g); // g = top of mut.perp. line 
	vmul(&y,e.y);
	vsub(d,y,h); // h = bot of mut.perp. line 
        if ((e.x>0.0 && e.x<1.0) && (e.y>0.0 && e.y<1.0)) { // in both lines 
		return 0;
	}
        if (e.y>0.0 && e.y<1.0) { // just in c-d 
		if(endOline(c,d,a)<endOline(c,d,b)) {
			vcopy(a,g); return 1;
		} else {
			vcopy(b,g); return 2;
		}
	}
        if (e.x>0.0 && e.x<1.0) { // just in a-b 
		if(endOline(a,b,c),endOline(a,b,d)) {
			vcopy(c,h); return 3;
		} else {
			vcopy(d,h); return 4;
		}
	}
	return -fsort4min(vdif(a,c),vdif(a,d),vdif(b,c),vdif(b,d));
}

float lineOline ( Vec a, Vec b, Vec c, Vec d, Vec *box ) {
// return overlap length for line segments a-b and c-d 
Vec x, y, z;
Mat M[1], W[1];
float lap, det, aa, ga, gb, hc, hd, r, s, t[4];
float ab,cd,ac,bd;
Vec e,f,g,h, pox[4],aox[4];
int i, key[4], ley[4];
int parr=0, perp=0;
	vcopy(a,aox+0); vcopy(b,aox+1); vcopy(c,aox+2); vcopy(d,aox+3);
	vcopy(a,box+0); vcopy(b,box+1); vcopy(c,box+2); vcopy(d,box+3);
	vsub(b,a,&x); vnorm(&x);
	vsub(d,c,&y); vnorm(&y);
	vprod(x,y,&z); vnorm(&z);
	if (fabs(vdot(x,y)) > 0.9999) { int in = 0;
        // parallel or anti
                if (in==0) in = dot2pair(a,b,c,&g,&h);
                if (in==0) in = dot2pair(a,b,d,&g,&h);
                if (in==0) in = dot2pair(c,d,a,&h,&g);
                if (in==0) in = dot2pair(c,d,b,&h,&g);
                if (in==0) return 0.0;
		parr=1;
        } else {
		VtoM(x,y,z,M);
		det = Mdet(M);
		Minv(M,W,det);
		vsub(d,a,&e);
		MmulV(W,e,&f);
		vmul(&x,f.x);
		vadd(a,x,&g); // g = top of mut.perp. line 
		vmul(&y,f.y);
		vsub(d,y,&h); // h = bot of mut.perp. line
	}
	ga = vdif(g,a); gb = vdif(g,b);
	hc = vdif(h,c); hd = vdif(h,d);
	vsub(a,g,&a); vsub(b,g,&b);
	vsub(c,h,&c); vsub(d,h,&d);
	if (vdot(a,c)==0.0) {
        // perpendicular
                if (fabs(ga-hc) > fabs(ga-hd)) hc = -hc;
		perp=1;
        } else {
                if (vdot(a,c)<0.0) hc = -hc;
        }
	if (vdot(a,b)<0.0) gb = -gb;
	if (vdot(a,d)<0.0) hd = -hd;
	// dists now relative to the contact normal gh: right+, left-
	t[0] = ga; t[1] = gb; t[2] = hc; t[3] = hd; 
	sort(0,t,0,key,4,1);
	for (i=0; i<4; i++) if (key[i]<2) ley[i] = 0; else ley[i] = 1;
	if (ley[0]==ley[1]) return 0.0;
	lap = fabs(t[key[1]]-t[key[2]]);
	if (box)
	{ int ke1 = key[1], ke2 = key[2],
	      le1 = ley[1], le2 = ley[2];
	  Vec tmp1, tmp2, tmp;
		if (ga > gb) vsub(box[0],box[1],&x);
			else vsub(box[1],box[0],&x);
		if (hc > hd) vsub(box[2],box[3],&y);
			else vsub(box[3],box[2],&y);
		vnorm(&x); vnorm(&y);
		vcopy(box[ke1],&tmp1);
		vcopy(box[ke2],&tmp2);
		vcopy(tmp1,box+0);
		if (le1==le2) { // contained 
			vcopy(tmp2,box+1);
			if (le2==0) {
				 vcopy(y,&tmp); vmul(&tmp,t[ke2]); vadd(h,tmp,box+3);
			} else { vcopy(x,&tmp); vmul(&tmp,t[ke2]); vadd(g,tmp,box+3); }
		} else { // staggered
			vcopy(tmp2,box+3);
			if (le2==0) {
				 vcopy(y,&tmp); vmul(&tmp,t[ke2]); vadd(h,tmp,box+1);
			} else { vcopy(x,&tmp); vmul(&tmp,t[ke2]); vadd(g,tmp,box+1); }
		}
		if (le1==0) {
			 vcopy(y,&tmp); vmul(&tmp,t[ke1]); vadd(h,tmp,box+2);
		} else { vcopy(x,&tmp); vmul(&tmp,t[ke1]); vadd(g,tmp,box+2); }
	}
	// restore ab cd line order in box (but box lines run parallel
	if (ley[1]) {
		pox[0] = box[2]; pox[1] = box[3];
		pox[2] = box[0]; pox[3] = box[1];
		for (i=0; i<4; i++) box[i] = pox[i];
	}
        ab = vdif(box[0],box[1]);
        cd = vdif(box[2],box[3]);
        if (fabs(ab-cd) > 0.0001) { float gh;
                ac = vdif(box[0],box[2]);
                bd = vdif(box[1],box[3]);
                printf("*NB* box has diff edges: ab=%f, cd=%f, ac=%f, bd=%f\n", ab,cd,ac,bd);
                Pv(aox[0]) NL Pv(aox[1]) NL Pv(aox[2]) NL Pv(aox[3]) NL Pv(h) NL Pv(g) NL
                ab = vdif(aox[0],aox[1]); cd = vdif(aox[2],aox[3]);
                ac = vdif(aox[0],aox[2]); bd = vdif(aox[1],aox[3]);
		gh = vdif(g,h);
                printf("     end-end distances : ab=%f, cd=%f, ac=%f, bd=%f, gh=%f\n", ab,cd,ac,bd,gh);
                printf("     end-box distances : ga=%f, gb=%f, hc=%f, hd=%f\n", ga,gb,hc,hd);
                ga = vdif(aox[0],g); gb = vdif(aox[1],g);
                hc = vdif(aox[0],h); hd = vdif(aox[3],h);
                printf("     end-box distances : ga=%f, gb=%f, hc=%f, hd=%f\n", ga,gb,hc,hd);
		if(perp) printf("perpendicular\n");
		if(parr) printf("parallel\n");
        }
	return lap;
}

float old_lineOline ( Vec a, Vec b, Vec c, Vec d, Vec *box ) {
// return overlap length for line segments a-b and c-d 
Vec x, y, z;
Mat M[1], W[1];
float lap, det, aa, ga, gb, hc, hd, r, s, t[4];
Vec e,f,g,h, pox[4];
int i, key[4], ley[4];
	vcopy(a,box+0); vcopy(b,box+1); vcopy(c,box+2); vcopy(d,box+3);
	vsub(b,a,&x); vnorm(&x);
	vsub(d,c,&y); vnorm(&y);
	vprod(x,y,&z); vnorm(&z);
	if (vdot(x,y) < 0.0001) { 
		vcopy(x,&e); vmul(&e,0.0001);
		vsub(c,e,&c); vadd(d,e,&d);
		vsub(d,c,&y); vnorm(&y);
	 }
	VtoM(x,y,z,M);
	det = Mdet(M);
	Minv(M,W,det);
	vsub(d,a,&e);
	MmulV(W,e,&f);
	vmul(&x,f.x);
	vadd(a,x,&g); // g = top of mut.perp. line
	vmul(&y,f.y);
	vsub(d,y,&h); // h = bot of mut.perp. line 
	ga = vdif(g,a); gb = vdif(g,b);
	hc = vdif(h,c); hd = vdif(h,d);
	vsub(a,g,&a); vsub(b,g,&b);
	vsub(c,h,&c); vsub(d,h,&d);
	aa = vsqr(a);
	r = vdot(a,c)/aa;
	s = vdot(a,d)/aa;
	if (vdot(a,b)<0.0) gb = -gb;
	if (r<0.0) hc = -hc;
	if (s<0.0) hd = -hd;
	t[0] = ga; t[1] = gb; t[2] = hc; t[3] = hd; 
	sort(0,t,0,key,4,1);
	for (i=0; i<4; i++) if (key[i]<2) ley[i] = 0; else ley[i] = 1;
	if (ley[0]==ley[1]) return 0.0;
	lap = fabs(t[key[1]]-t[key[2]]);
	if (box)
	{ int ke1 = key[1], ke2 = key[2],
	      le1 = ley[1], le2 = ley[2];
	  Vec tmp1, tmp2, tmp;
		if (ga > gb) vsub(box[0],box[1],&x);
			else vsub(box[1],box[0],&x);
		if (hc > hd) vsub(box[2],box[3],&y);
			else vsub(box[3],box[2],&y);
		vnorm(&x); vnorm(&y);
		vcopy(box[ke1],&tmp1);
		vcopy(box[ke2],&tmp2);
		vcopy(tmp1,box+0);
		if (le1==le2) { // contained 
			vcopy(tmp2,box+1);
			if (le2==0) {
				 vcopy(y,&tmp); vmul(&tmp,t[ke2]); vadd(h,tmp,box+3);
			} else { vcopy(x,&tmp); vmul(&tmp,t[ke2]); vadd(g,tmp,box+3); }
		} else { // staggered
			vcopy(tmp2,box+3);
			if (le2==0) {
				 vcopy(y,&tmp); vmul(&tmp,t[ke2]); vadd(h,tmp,box+1);
			} else { vcopy(x,&tmp); vmul(&tmp,t[ke2]); vadd(g,tmp,box+1); }
		}
		if (le1==0) {
			 vcopy(y,&tmp); vmul(&tmp,t[ke1]); vadd(h,tmp,box+2);
		} else { vcopy(x,&tmp); vmul(&tmp,t[ke1]); vadd(g,tmp,box+2); }
	}
	// restore ab cd line order in box (but box lines run parallel
	if (ley[1]) {
		pox[0] = box[2]; pox[1] = box[3];
		pox[2] = box[0]; pox[3] = box[1];
		for (i=0; i<4; i++) box[i] = pox[i];
	}
	return lap;
}

int line2dot ( Vec a, Vec b, Vec c, float s ) {
// TRUE if c lies over line segment a-b closer than s
Vec p, q;
float d;
	vsub(b,a,&p);
	vsub(c,a,&q);
	d = vdot(p,q)/vsqr(p);
	if (d < 0.0 || d > 1.0) return 0;
	vmul(&p,d);
	vsub(q,p,&q);
	d = vsqr(q);
	if (d>s*s) return 0;
	return 1;
}

float endOline ( Vec a, Vec b, Vec c ) {
// returns distance of c to an extended line a-b 
Vec p, q;
float d;
	vsub(b,a,&p);
	vsub(c,a,&q);
	d = vdot(p,q)/vsqr(p);
	vmul(&p,d);
	vsub(q,p,&q);
	return vmod(q);
}

float dotOline ( Vec a, Vec b, Vec c, Vec *e ) {
// returns distance of c to an extended line a-b and image of c on a-b in e 
Vec p, q;
float d;
	vsub(b,a,&p);
	vsub(c,a,&q);
	d = vdot(p,q)/vsqr(p);
	vmul(&p,d);
	vsub(q,p,&q);
	d = vsqr(q);
	vadd(a,p,&p);
	if (e) { e->x = p.x; e->y = p.y; e->z = p.z; }
	return sqrt(d);
}
*/
