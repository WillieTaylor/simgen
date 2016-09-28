float randf();

class Vec {
public:
	float x, y, z;
	// constructors
	Vec(): x(0), y(0), z(0) {}
	Vec( float a,float b,float c ): x(a), y(b), z(c) {}
	Vec( float s ): x(s*2.0*(randf()-0.5)), y(s*2.0*(randf()-0.5)), z(s*2.0*(randf()-0.5)) {}
	// operators
	Vec operator + ( const Vec &v ) const {		// add two vectors
		return Vec( x+v.x, y+v.y, z+v.z );
	}
	Vec operator + ( const float s ) const {	// add a scalar to each component
		return Vec( x+s, y+s, z+s );
	}
	Vec operator - ( const Vec &v ) const {		// subtract two vectors
		return Vec( x-v.x, y-v.y, z-v.z );
	}
	Vec operator - ( const float s ) const {	// subtract a scalar from each component
		return Vec( x-s, y-s, z-s );
	}
	Vec operator - ()  {				// negate
		return Vec(-x, -y, -z);
	}
	Vec operator * ( const float s ) const {	// multiply each component by a scalar
		return Vec( x*s, y*s, z*s );
	}
	void operator *= ( const float s ) {	// multiply each component by a scalar
		x*=s; y*=s; z*=s;
	}
	float operator * ( const Vec &v ) const {	// scalar product
		return x*v.x + y*v.y + z*v.z;
	}
	Vec operator ^ ( const Vec &v ) const {		// vector product
		return Vec( y*v.z - v.y*z, z*v.x - v.z*x, x*v.y - v.x*y);
	}
	Vec operator *  ( const Mat &M );		// multiply a Vec and a Mat
	Vec operator *= ( const Mat &M );		// multiply a Vec and a Mat and set
		// defined after Mat class declared
	Vec operator / ( const float s ) const {	// divide each component by a scalar
		return Vec( x/s, y/s, z/s );
	}
	void operator /= ( const float s ) {	// multiply each component by a scalar
		x/=s; y/=s; z/=s;
	}
	void operator += ( const Vec &v ) {		// sum a vector into the current vector
		x += v.x; y += v.y; z += v.z;
	}
	void operator -= ( const Vec &v ) {		// subtract a vector from the current vector
		x = x-v.x; y = y-v.y; z = z-v.z;
	}
	float operator | ( const Vec &v ) const {	// distance between two vectors
		float a = x-v.x, b = y-v.y, c = z-v.z;
		return sqrtf(a*a + b*b + c*c);
	}	
	float operator || ( const Vec &v ) const {	// squared distance between two vectors
		float a = x-v.x, b = y-v.y, c = z-v.z;
		return a*a + b*b + c*c;
	}	
	Vec operator & ( const Vec &v ) const {		// average of two vectors
		return Vec( (x+v.x)*0.5, (y+v.y)*0.5, (z+v.z)*0.5 );
	}
	Vec operator &= ( const Vec &v ) {		// average given vector with current
		x = (x+v.x)*0.5; y = (y+v.y)*0.5; z = (z+v.z)*0.5;
		return *this;
	}	
	// inline functions
	inline bool iszero () { if (sqrt(x*x+y*y+z*z) < NOISE) return 1; else return 0; }
	inline bool isunit () { if (sqrt(x*x+y*y+z*z-1.0) < NOISE) return 1; else return 0; }
	inline void zero () { x = y = z = 0.0; }	
	inline void init () { x = y = z = 0.0; }	
	inline float len () const {			// vector length
		return sqrtf(x*x + y*y + z*z);
	}
	inline float mod () const {			// vector length
		return sqrtf(x*x + y*y + z*z);
	}
	inline float sqr () const {			// vector length squared
		return x*x + y*y + z*z;
	}	
	inline float ddist ( const Vec &v ) const {	// distance between two vectors
		float a = x-v.x, b = y-v.y, c = z-v.z;
		return a*a + b*b + c*c;
	}	
	inline float dist ( const Vec &v ) const {	// squared distance between two vectors
		return sqrtf(ddist(v));
	}	
	inline void set_ave ( const Vec &v ) {		// average of current and given vectors
		x = (x+v.x)*0.5; y = (y+v.y)*0.5; z = (z+v.z)*0.5;
	}
	inline void setVec () {				// normalise the current vector
		float s = sqrtf(x*x + y*y + z*z);
		x/=s; y/=s; z/=s;
	}
	inline void setVec ( const float len ) {	// normalise the current vector to given length
		float s = len/sqrtf(x*x + y*y + z*z);
		x*=s; y*=s; z*=s;
	}
	inline void setVec ( const Vec &v ) {		// set to normalised given vector (unchanged)
		float s = v.len();
		x=v.x/s; y=v.y/s; z=v.z/s;
	}
	inline Vec getVec () const {			// return the normalised current vector
		float s = sqrtf(x*x + y*y + z*z);
		return Vec( x/s, y/s, z/s );
	}
	inline Vec getVec ( const float len ) const {	// return the current vector with given length
		float s = len/sqrtf(x*x + y*y + z*z);
		return Vec( x*s, y*s, z*s );
	}
	inline Vec norm () const {			// return the normalised current vector
		float s = sqrtf(x*x + y*y + z*z);
		return Vec( x/s, y/s, z/s );
	}
	// a.setVec(); sets a to be unit length
	// a.setVec(b); sets a to the unit vector along b (with b unchanged)
	// b = a.getVec(); sets b to the unit vector along a (with a unchanged)
	// b = a.norm(); same as above
	//
	// functions
	void set_disp ();			// add a random displacement between +/-1 to each component
	void set_disp ( const float s );	// add a random displacement between +/-<s> (was vradd())
	void set_rand ();			// set a random vector between +/-1 to each component
	void set_rand ( const float s );	// set a random vector between +/-<s> (was vrand())
	bool	vec_in_seg  ( const Seg& ) const; 	// TRUE if c lies over the line segment
	bool	vec_in_seg  ( const Vec&, const Vec&) const; 	// dito
	float	vec_to_line ( const Vec& ) const;	// returns the distance of Vec to extended 0-B line segment
	float	vec_to_line ( const Seg& ) const;	// returns the distance of Vec to extended A-B line segment
	float	vec_to_line ( const Vec&, const Vec& ) const;	// dito
	Vec	vec_on_line ( const Seg& ) const;	// returns the image of Vec on extended line segment
	Vec	vec_on_line ( const Vec&, const Vec& ) const;	// dito
	Vec  get_rot ( const Seg&, const float ) const;	// rotate about line (Seg) by angle
	Vec  get_rot ( const Vec&, const Vec&, const float ) const; // dito
	void set_rot ( const Seg&, const float );	// rotate about line (Seg) by angle
	void set_rot ( const Vec&, const Vec&, const float ); //dito with 2 vectors
	void set_rot ( const Vec&, const float ); //dito with origin
	void set_rot ( const float ); //dito with origin and random
	Vec  get_frac ( const Vec&, const Mat& );	// returns the components of Vec in basis-set Mat
	Vec  get_frac ( const Mat& );		// returns the components of <this> Vec in basis-set Mat
	void set_frac ( const Mat& );		// sets the components of <this> Vec in basis-set Mat
	void set_frac ( const Vec&, const Mat& );	// dito with given Vec
};
