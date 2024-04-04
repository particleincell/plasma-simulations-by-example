/*defines the simulation domain*/
#include <random>
#include <math.h>
#include <iostream>
#include "World.h"
#include "Field.h"
#include "Species.h"
	
//make an instance of the Rnd class
Rnd rnd;

using namespace std;

/*constructor*/
World::World(int ni, int nj):
	ni{ni}, nj{nj},	nn{ni,nj},
	phi(nn),rho(nn),node_vol(nn),
	ef(nn), object_id(nn)	{
		time_start =  chrono::high_resolution_clock::now();	//save starting time point
	}

/*sets domain bounding box and computes mesh spacing*/
void World::setExtents(double2 _x0, double2 _xm) {
	/*set origin and the opposite corner*/
	x0 = _x0;
	xm = _xm;

	/*compute spacing by dividing length by the number of cells*/
	for (int i=0;i<2;i++)
		dh[i] = (xm(i)-x0(i))/(nn(i)-1);

	//compute centroid
	xc = 0.5*(x0+xm);

	/*recompute node volumes*/
	computeNodeVolumes();
}

/*returns elapsed wall time in seconds*/
double World::getWallTime() {
  auto time_now = chrono::high_resolution_clock::now();
  chrono::duration<double> time_delta = time_now-time_start;
  return time_delta.count();
}

/*computes charge density from rho = sum(charge*den)*/
void World::computeChargeDensity(vector<Species> &species)
{
	rho = 0;
	for (Species &sp:species)
	{
		if (sp.charge==0) continue;	//don't bother with neutrals
		rho += sp.charge*sp.den;
	}
}	

/*computes node volumes on axisymmetric domain*/
void World::computeNodeVolumes() {
	for (int i=0;i<ni;i++)
		for (int j=0;j<nj;j++)
			{
				double r2 = getR(min(j+0.5,nj-1.0));
				double r1 = getR(max(j-0.5,0.0));
				double dz = dh[0];
				if (i==0 || i==ni-1) dz*=0.5;
				double V= dz*Const::PI*(r2*r2-r1*r1);
				if (j==0) V*=2.0;
				node_vol[i][j] = V;
			}
}

/* computes total potential energy from 0.5*eps0*sum(E^2)*/
double World::getPE() {
	double pe = 0;
	for (int i=0;i<ni;i++)
		for (int j=0;j<nj;j++)
			{
				double3 efn = ef[i][j];	//ef at this node
				double ef2 = efn[0]*efn[0]+efn[1]*efn[1];

				pe += ef2*node_vol[i][j];
			}
	return 0.5*Const::EPS_0*pe;
}

/*sugarcubes a circle centered at (x0,y0)*/
void World::addCircle(double2 x0, double radius, double phi_circle)
{
    /*save circle parameters*/
    circle_x0 = x0;
    circle_rad2 = radius*radius;

    for (int i=0;i<ni;i++)
        for (int j=0;j<nj;j++)
        {
            /*compute node position*/
            double3 x = pos(i,j);
            if (inCircle(x))
            {
                object_id[i][j] = 1;
                phi[i][j] = phi_circle;
            }
        }
}

/*marks k=0 plane as 0V Dirichlet boundary*/
void World::addInlet() {
	for (int j=0;j<nj;j++)
	{
		object_id[0][j] = 2;
		phi[0][j] = 0;
	}
}
	
/*returns true if point x is inside or on the circle*/
bool World::inCircle(double3 x)
{
	double2 x2 {x(0),x(1)};
	double2 r = x2-circle_x0;	//ray to x

    double r_mag2 = (r[0]*r[0] + r[1]*r[1]);
    if (r_mag2<=circle_rad2) return true;
    return false;
}

/*checks for steady state by comparing change in mass, momentum, and energy*/
bool World::steadyState(vector<Species> &species) {
	// do not do anything if already at steady state
	if (steady_state) return true;

	double tot_mass = 0;
	double tot_mom = 0;
	double tot_en = getPE();
	for (Species &sp:species)
	{
		tot_mass += sp.getRealCount();	//number of real molecules
		tot_mom += mag(sp.getMomentum());	//momentum magnitude
		tot_en += sp.getKE();		//add kinetic energy
	}

	/*compute new values to last*/
	const double tol = 1e-3;
	if (abs((tot_mass-last_mass)/tot_mass)<tol &&
		abs((tot_mom-last_mom)/tot_mom)<tol &&
		abs((tot_en-last_en)/tot_en)<tol) {
		steady_state = true;
		cout<<"Steady state reached at time step "<<ts<<endl;
	}

	/*update prior values*/
	last_mass = tot_mass;
	last_mom = tot_mom;
	last_en = tot_en;
	return steady_state;
}

