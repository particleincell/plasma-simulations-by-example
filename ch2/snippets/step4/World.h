#ifndef _WORLD_H
#define _WORLD_H

#include <random>
#include <vector>
#include "Field.h"

class Species;

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

//World.h
class World {
public:
  World (int ni, int nj, int nk);  //constructor
  
  //sets the mesh span, also recomputes cell spacing
  void setExtents(const double3 x0, const double3 xm);

  double3 getX0() const {return double3(x0);}
  double3 getXm() const {return double3(xm);}
  double3 getXc() const {return double3(xc);}
  double3 getDh() const {return double3(dh);}

  /*converts physical position to logical coordinate*/
  double3 XtoL(double3 x) const {
	 double3 lc;
	 lc[0] = (x[0]-x0[0])/dh[0];
	 lc[1] = (x[1]-x0[1])/dh[1];
	 lc[2] = (x[2]-x0[2])/dh[2];
	 return lc;
  }

  /*computes rho from the sum of species charge and density*/    
  void computeChargeDensity(std::vector<Species> &species);

  const int nn[3];	 	//number of nodes  
  const int ni,nj,nk;	//number of nodes in individual variables
  
  Field phi;		//potential
  Field rho;		//charge density
  Field node_vol;	//node volumes
  Field3 ef;		//electric field

protected:
  double x0[3];	 	//mesh origin
  double dh[3];	 	//cell spacing
  double xm[3];		//mesh max bound
  double xc[3];		//domain centroid

  void computeNodeVolumes();
};

#endif
