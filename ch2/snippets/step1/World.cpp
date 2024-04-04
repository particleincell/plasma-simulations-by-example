//World.cpp
#include "World.h"

//constructor
World::World (int ni, int nj, int nk) : 
	ni{ni}, nj{nj}, nk{nk},	nn{ni,nj,nk} {}
	
//sets mesh extents and computes cell spacing
void World::setExtents(double x1, double y1, double z1, 
					double x2, double y2, double z2) 
{
	/*set origin*/
	x0[0] = x1; x0[1] = y1; x0[2] = z1;

	/*opposite corner*/
	xm[0] = x2; xm[1] = y2; xm[2] = z2;

	/*compute spacing by dividing length by number of cells*/
	for (int i=0;i<3;i++)
		dh[i] = (xm[i]-x0[i])/(nn[i]-1);

	for (int i=0;i<3;i++)			//compute centroid
		xc[i] = 0.5*(x0[i]+xm[i]);

	/*recompute node volumes*/
			computeNodeVolumes();
}

void World::computeNodeVolumes() 
{
	/*to be implemented*/
}
