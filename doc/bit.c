float fixgeom ( Cell *cell, float fix )
	:
// quality checks
	d = vdif(b1->xyz,a1->xyz);
	dlo = c0->prox.y*(1.0-good);
	dhi = c0->prox.y*(1.0+good);
	if (d<dlo || d>dhi) return 99.9; // proximal distance too poor
	d = vdif(b2->xyz,a2->xyz);
	dlo = c0->dist.y*(1.0-good);
	dhi = c0->dist.y*(1.0+good);
	if (d<dlo || d>dhi) return 99.9; // distal distance too poor
	d = angle(b2->xyz,c0->xyz,a2->xyz)*180/PI;
	if (d < 30.0 || d > 150.0) return 0.0; // base triangle too thin
	t1 = torsion(b2->xyz,b1->xyz,c0->xyz,a1->xyz);
	if (t1>999.0) { Pt(Bad t1 in) Pi(cell->uid) NL
		Pv(b2->xyz) NL Pv(b1->xyz) NL Pv(c0->xyz) NL Pv(a1->xyz) NL exit(1);
	}
	th = angle(b1->xyz, c0->xyz, a1->xyz);
	t2 = torsion(a2->xyz,a1->xyz,c0->xyz,b1->xyz);
	if (t2>999.0) { Pt(Bad t2 in) Pi(cell->uid) NL
		Pv(a2->xyz) NL Pv(a1->xyz) NL Pv(c0->xyz) NL Pv(b1->xyz) NL exit(1);
	}
	dt1old = angdif(t1,tau1);
	dthold = angdif(th,theta);
	dt2old = angdif(t2,tau2);
	if (dt1old+dthold+dt2old > 3.0) return 99.9; // too far to fix
// configuration close enough to fix
	xyz[0] = b2->xyz; xyz[1] = b1->xyz;
	xyz[2] = c0->xyz;
	xyz[3] = a1->xyz; xyz[4] = a2->xyz;
	// dial-up correct torsion angles
	t1 = torsion(xyz[3],xyz[2],xyz[1],xyz[0]);
	xyz[0].get_rot(xyz[2],xyz[1],tau1-t1);
	t2 = torsion(xyz[1],xyz[2],xyz[3],xyz[4]);
	xyz[4].get_rot(xyz[2],xyz[3],tau2-t2);
	t1 = torsion(xyz[0], xyz[1], xyz[2], xyz[3]);
	th = angle(xyz[1], xyz[2], xyz[3]);
	t2 = torsion(xyz[1], xyz[2], xyz[3], xyz[4]);
	dt1new = angdif(t1,tau1);
	dthnew = angdif(th,theta);
	dt2new = angdif(t2,tau2);
	// fit positions xyz[0,2,4] to b2,c0,a2
	x = xyz[4]- xyz[0];
	mid = xyz[4] & xyz[0];
	y = xyz[2] - mid;
	z = x^y;
	mat = Mat(x.norm(), y.norm(), z.norm()); // new unit basis vectors in mat
	wat = mat.get_inv();
	mid = (xyz[0] + xyz[2] + xyz[4])/3.0;	// new CoG
	for (i=0; i<5; i++) {  // get xyz from centre in terms of basis set coeficients (abc)
		xyz[i] -= mid;
		abc[i] = wat * xyz[i];
	}
	x = a2->xyz - b2->xyz;
	mid = a2->xyz & b2->xyz;
	y = c0->xyz - mid;
	z = x^y;
	x.setVec(); y.setVec(); z.setVec();		// old basis vectors
	mid = (b2->xyz + c0->xyz + a2->xyz)/3.0;	// old CoG
	for (i=0; i<5; i++) {	// reconstruct new positions with old basis vectors
		xyz[i] = mid;	// start at centre and sum basis vector components
		xyz[i] += x * abc[i].x;
		xyz[i] += y * abc[i].y;
		xyz[i] += z * abc[i].z;
	}

	t1 = torsion(xyz[0], xyz[1], xyz[2], xyz[3]);
	th = angle(xyz[1], xyz[2], xyz[3]);
	t2 = torsion(xyz[1], xyz[2], xyz[3], xyz[4]);
	dt1new = angdif(t1,tau1);
	dthnew = angdif(th,theta);
	dt2new = angdif(t2,tau2);
	if (dt1new>dt1old || dthnew>dthold || dt2new>dt2old) return 0.0;
// shift towards new positions (by fix) if all angles better
	shift = (xyz[0]-b2->xyz)*fix; moveCell(b2,shift);
	shift = (xyz[1]-b1->xyz)*fix; moveCell(b1,shift);
	shift = (xyz[2]-c0->xyz)*fix; moveCell(c0,shift);
	shift = (xyz[3]-a1->xyz)*fix; moveCell(a1,shift);
	shift = (xyz[4]-a2->xyz)*fix; moveCell(a2,shift);
	t1 = torsion(b2->xyz,b1->xyz,c0->xyz,a1->xyz);
	th = angle(b1->xyz, c0->xyz, a1->xyz);
	t2 = torsion(a2->xyz,a1->xyz,c0->xyz,b1->xyz);
	t1 = angdif(t1,tau1);
	th = angdif(th,theta);
	t2 = angdif(t2,tau2);
	tt = t1 + th + t2;
	return tt;
}
