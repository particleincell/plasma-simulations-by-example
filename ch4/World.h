#ifndef _WORLD_H
#define _WORLD_H

#include <vector>
#include <random>
#include <chrono>
#include "Field.h"

class Species;

/*define constants*/
namespace Const
{
	const double EPS_0 = 8.85418782e-12;  	// C/(V*m), vacuum permittivity
	const double QE = 1.602176565e-19;		// C, electron charge
	const double AMU = 1.660538921e-27;		// kg, atomic mass unit
	const double ME = 9.10938215e-31;		// kg, electron mass
	const double K = 1.380648e-23;			// J/K, Boltzmann constant
	const double PI = 3.141592653;			// pi
	const double EvToK = QE/K;				// 1eV in K ~ 11604
}

/*object for sampling random numbers*/
class Rnd {
public:
	//constructor: set initial random seed and distribution limits
	Rnd(): mt_gen{std::random_device()()}, rnd_dist{0,1.0} {}
	double operator() () {return rnd_dist(mt_gen);}

protected:
	std::mt19937 mt_gen;	    //random number generator
	std::uniform_real_distribution<double> rnd_dist;  //uniform distribution
};

extern Rnd rnd;		//tell the compiler that an object of type Rnd called rnd is defined somewhere

/*defines the computational domain*/
class World
{
public:	
	/*constructor, allocates memory*/
	World(int ni, int nj, int nk);

	/*functions to set mesh origin and spacing*/
	void setExtents(const double3 x0, const double3 xm);
	
	double3 getX0() const {return double3(x0);}
	double3 getXm() const {return double3(xm);}
	double3 getXc() const {return double3(xc);}
	double3 getDh() const {return double3(dh);}
	double getCellVolume() const {return dh(0)*dh(1)*dh(2);} //cell volume

	/*functions for accessing time information*/
	int getTs() const {return ts;}
	double getTime() const {return time;}
	double getWallTime();  /*returns wall time in seconds*/
	double getDt() const {return dt;}
	bool isLastTimeStep() const {return ts==num_ts-1;}

	bool inBounds(double3 pos) {
		for (int i=0;i<3;i++)
			if (pos[i]<x0[i] || pos[i]>=xm[i]) return false;
		return true;
	}

	/*sets time step and number of time steps*/
	void setTime(double dt, int num_ts) {this->dt=dt;this->num_ts=num_ts;}
	
	/*advances to the next time step, returns true as long as more time steps remain*/
	bool advanceTime() {time+=dt;ts++;return ts<=num_ts;}

	/*checks and sets a steady state flag*/
	bool steadyState(std::vector<Species> &species);

	/*returns steady state flag*/
	bool isSteadyState() {return steady_state;}

	/*converts physical position to logical coordinate*/
	double3 XtoL(double3 x) const {
	  	double3 lc;
		lc[0] = (x[0]-x0(0))/dh(0);
		lc[1] = (x[1]-x0(1))/dh(1);
		lc[2] = (x[2]-x0(2))/dh(2);
		return lc;
	}

	/*returns cell index for a given point*/
	int3 XtoIJK(double3 x) const {
		double3 lc = XtoL(x);
		int3 ijk {(int)lc[0],(int)lc[1],(int)lc[2]};
		return ijk;
		}

	/*returns cell index for a position x*/
	int XtoC(double3 x) const {
		int3 i = XtoIJK(x);
		return i[2]*(nj-1)*(ni-1)+i[1]*(ni-1)+i[0];
	}


	/*converts logical coordinate to physical position*/
	double3 pos(double3 lc)
	{
		double3 x = x0+dh*lc;
		return x;
	}

	/*another form that takes 3 ints as inputs*/
	double3 pos(int i, int j, int k) {
		double3 x{(double)i,(double)j,(double)k};
		return pos(x);
	}

	int U(int i,int j, int k) {return object_id.U(i,j,k);}

	/*computes charge density from rho = sum(charge*den)*/
	void computeChargeDensity(std::vector<Species> &species);

	/*returns the system potential energy*/
	double getPE();

	/*sugarcubes a sphere centered at (x0,y0,z0)*/
	void addSphere(double3 x0, double radius, double phi_sphere);
	
	/*return true if point x is inside or on the sphere*/
	bool inSphere(double3 x);

	/*marks the k=0 face as Dirichlet*/
	void addInlet();

	/*returns the parametric position for the intersection point*/
	double lineSphereIntersect(const double3 &x1, const double3 &x2);

	/*returns random diffuse velocity sampled at a given surface point*/
	double3 sphereDiffuseVector(const double3 &x);

	//mesh geometry
	const int ni,nj,nk;	//number of nodes
	const int3 nn;	//another way to access node counts

	Field phi;			//potential
	Field rho;			//charge density
	Field node_vol;		//node volumes
	Field3 ef;			//electric field components
	FieldI object_id; 	//object id flag to flag fixed nodes

protected:
	double3 x0;	//mesh origin
	double3 dh;	//cell spacing

	double3 xm;		//origin-diagonally opposite corner (max bound)
	double3 xc;		//domain centroid
	
	double dt = 0;		//time step length
	double time = 0;	//physical time
	int ts = -1;		//current time step
	int num_ts = 0;		//number of time steps

	/*sphere data*/
	double3 sphere_x0 {0,0,0};	//sphere centroid
	double sphere_rad2 = 0;		//sphere radius squared

	std::chrono::time_point<std::chrono::high_resolution_clock> time_start;	//time at simulation start

	bool steady_state = false;	//set to true once steady state is reached
	double last_mass = 0;	//mass at the prior time step
	double last_mom = 0;	//momentum at the prior time step
	double last_en = 0;		//energy at the prior time step

	void computeNodeVolumes();
};

#endif
