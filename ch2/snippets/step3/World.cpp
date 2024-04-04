//World.cpp
#include "World.h"
#include "Species.h"

Rnd rnd;			//create an instance of a Rnd object

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

	/*recompute node volumes*/
	computeNodeVolumes(); 
}

/*computes node volumes, dx*dy*dz on internal nodes and fractional
 * values on domain boundary faces*/
void World::computeNodeVolumes() {
	for (int i=0;i<ni;i++)
		for (int j=0;j<nj;j++)
			for (int k=0;k<nk;k++)
			{
				double V = dh[0]*dh[1]*dh[2];	//default volume
				if (i==0 || i==ni-1) V*=0.5;	//reduce by two for each boundary index
				if (j==0 || j==nj-1) V*=0.5;
				if (k==0 || k==nk-1) V*=0.5;
				node_vol[i][j][k] = V;
			}
}
 
/*computes charge density from rho = sum(charge*den)*/
void World::computeChargeDensity(std::vector<Species> &species)
{
	rho = 0;
	for (Species &sp:species)
	{
		if (sp.charge==0) continue;	//don't bother with neutrals
		rho += sp.charge*sp.den;
	}
}	
