class Mat {
public:
	Vec A, B, C;
	// constructors
	Mat(): A(Vec(1,0,0)), B(Vec(0,1,0)), C(Vec(0,0,1)) {}
	Mat( const Vec &a, const Vec &b, const Vec &c ): A(a), B(b), C(c) {}
	Mat( const Mat &M ): A(M.A), B(M.B), C(M.C) {}
	// operators
	Mat operator *  ( const Mat &Q );		// multiply two matrices
	Vec operator *  ( const Vec &v );		// multiply a Mat and a Vec
	Mat operator *= ( const Mat &Q );		// multiply and set matrices
	//
	// inline functions
	inline Mat print () const {
        	printf("\n");
		printf("%9.6f %9.6f %9.6f\n", A.x,A.y,A.z);
		printf("%9.6f %9.6f %9.6f\n", B.x,B.y,B.z);
		printf("%9.6f %9.6f %9.6f\n", C.x,C.y,C.z);
        	printf("\n");
	}
	inline void zero () {				// set all components to sero
		A.x=A.y=A.z = B.x=B.y=B.z = C.x=C.y=C.z = 0.0;
	}
	inline void ident () {				// set the identity matrix
		A.y=A.z = B.x=B.z = C.x=C.y = 0.0;
		A.x = B.y = C.z = 1.0;
	}
	inline float det () const {			// determinant of current matrix
		return A.x * (B.y*C.z - B.z*C.y)
        	     - B.x * (A.y*C.z - A.z*C.y)
        	     + C.x * (A.y*B.z - A.z*B.y);
	}
	inline Mat get_trans () const {			// set to the transpose of the current matrix
		Mat T, M = *this;
		T.A.x = M.A.x; T.A.y = M.B.x; T.A.z = M.C.x;
		T.B.x = M.A.y; T.B.y = M.B.y; T.B.z = M.C.y;
		T.C.x = M.A.z; T.C.y = M.B.z; T.C.z = M.C.z;
		return Mat(T);
	}
	inline void set_trans () {			// set the current matrix to its transpose
		Mat M = *this; // not pointer as this is overwritten
		A.x = M.A.x; A.y = M.B.x; A.z = M.C.x;
		B.x = M.A.y; B.y = M.B.y; B.z = M.C.y;
		C.x = M.A.z; C.y = M.B.z; C.z = M.C.z;
	}
	inline void set_trans ( const Mat &M ) {	// set the current matrix to M transpose
		A.x = M.A.x; A.y = M.B.x; A.z = M.C.x;
		B.x = M.A.y; B.y = M.B.y; B.z = M.C.y;
		C.x = M.A.z; C.y = M.B.z; C.z = M.C.z;
	}
	inline bool iszero () {
		if ( A.sqr() + B.sqr() + C.sqr() < NOISE ) return 1; else return 0;
	}
	inline bool isident () { if (A.x*A.x+B.y*B.y+C.z*C.z-3.0 < NOISE
		&& A.y*A.y+A.z*A.z + B.x*B.x+B.z*B.z + C.x*C.x+C.y*C.y < NOISE ) return 1; else return 0;
	}
	// functions
	void set_inv ( const float );		// set the current matrix to its inverse (given det)
	void set_inv ();			// set the current matrix to its inverse (find det)
	Mat get_inv ( const float ) const;	// return the inverse of the current matrix (given det)
	Mat get_inv () const;			// return the inverse of the current matrix (find det)
	Mat inv ( const float ) const;		// same
	Mat inv () const;			// same
	void set_rot ( const char axis, const float r ); // rotates the current matrix about axis
};
