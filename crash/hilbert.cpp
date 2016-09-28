	struct Turtle {
	    void Move(float dx, float dy, bool changeFace = false) {
	        if (changeFace) ++Face;
	        switch (Face) {
	            case 0: P += Vector3(dx, dy, 0); break;
	            case 1: P += Vector3(0, dy, dx); break;
	            case 2: P += Vector3(-dy, 0, dx); break;
	            case 3: P += Vector3(0, -dy, dx); break;
	            case 4: P += Vector3(dy, 0, dx); break;
	            case 5: P += Vector3(dy, dx, 0); break;
	        }
	        HilbertPath.push_back(P);
	    }
	    Point3 P;
	    int Face;
	};
	 
	static void HilbertU(int level)
	{
	    if (level == 0) return;
	    HilbertD(level-1);      HilbertTurtle.Move(0, -dist);
	    HilbertU(level-1);      HilbertTurtle.Move(dist, 0);
	    HilbertU(level-1);      HilbertTurtle.Move(0, dist);
	    HilbertC(level-1);
	}
	  
	static void HilbertD(int level)
	{
	    if (level == 0) return;
	    HilbertU(level-1);      HilbertTurtle.Move(dist, 0);
	    HilbertD(level-1);      HilbertTurtle.Move(0, -dist);
	    HilbertD(level-1);      HilbertTurtle.Move(-dist, 0);
	    HilbertA(level-1);
	}
	  
	static void HilbertC(int level)
	{
	    if (level == 0) return;
	    HilbertA(level-1);      HilbertTurtle.Move(-dist, 0);
	    HilbertC(level-1);      HilbertTurtle.Move(0, dist);
	    HilbertC(level-1);      HilbertTurtle.Move(dist, 0);
	    HilbertU(level-1);
	}
	  
	static void HilbertA(int level)
	{
	    if (level == 0) return;
	    HilbertC(level-1);      HilbertTurtle.Move(0, dist);
	    HilbertA(level-1);      HilbertTurtle.Move(-dist, 0);
	    HilbertA(level-1);      HilbertTurtle.Move(0, -dist);
	    HilbertD(level-1);
	}
	 
	void CreateHilbertCube(int lod)
	{
	    HilbertU(lod); HilbertTurtle.Move(dist, 0, true);
	    HilbertU(lod); HilbertTurtle.Move(0, dist, true);
	    HilbertC(lod); HilbertTurtle.Move(0, dist, true);
	    HilbertC(lod); HilbertTurtle.Move(0, dist, true);
	    HilbertC(lod); HilbertTurtle.Move(dist, 0, true);
	    HilbertD(lod);
	}
