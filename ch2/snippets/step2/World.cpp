//World.cpp
#include "World.h"

//constructor
World::World (int ni, int nj, int nk) : 
	ni{ni}, nj{nj}, nk{nk},	nn{ni,nj,nk},
	phi(ni,nj,nk),rho(ni,nj,nk),node_vol(ni,nj,nk),
	ef(ni,nj,nk)	{ }
	
//sets mesh extents and computes cell spacing
void World::setExtents(double3 _x0, double3 _xm) {
	/*set origin*/
	for (int i=0;i<3;i++)
		x0[i] = _x0[i];

	/*opposite corner*/
	for (int i=0;i<3;i++)
		xm[i] = _xm[i];

	/*compute spacing by dividing length by number of cells*/
	for (int i=0;i<3;i++)
		dh[i] = (xm[i]-x0[i])/(nn[i]-1);

	for (int i=0;i<3;i++)			//compute centroid
		xc[i] = 0.5*(x0[i]+xm[i]);
}

