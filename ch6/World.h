#ifndef _WORLD_H
#define _WORLD_H

#include <vector>
#include <random>
#include <chrono>
#include "Field.h"
#include "VolumeMesh.h"

class Species;
class Particle;

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
	World(VolumeMesh &vm);  /*constructor*/

	/*functions for accessing time information*/
	int getTs() const {return ts;}
	double getTime() const {return time;}
	double getWallTime();  /*returns wall time in seconds*/
	double getDt() const {return dt;}
	bool isLastTimeStep() const {return ts==num_ts-1;}

	/*pass through wrapper*/
	LCord XtoL(const double3 &pos, int e, bool search=true) {return vm.XtoL(pos,e,search);}

	/*brute force search*/
	LCord XtoLbrute(const double3 &post);

	/*sets time step and number of time steps*/
	void setTime(double dt, int num_ts) {this->dt=dt;this->num_ts=num_ts;}
	
	/*advances to the next time step, returns true as long as more time steps remain*/
	bool advanceTime() {time+=dt;ts++;return ts<=num_ts;}

	/*checks and sets a steady state flag*/
	bool steadyState(std::vector<Species> &species);

	/*computes charge density from rho = sum(charge*den)*/
	void computeChargeDensity(std::vector<Species> &species);

	/*returns the system potential energy*/
	double getPE();

	//mesh geometry
	VolumeMesh &vm;		//volume mesh
	Field node_vol;		//node volumes
	Field phi;			//potential
	Field rho;			//charge density
	Field3 ef;			//electric field components

protected:
	double dt = 0;		//time step length
	double time = 0;	//physical time
	int ts = -1;		//current time step
	int num_ts = 0;		//number of time steps

	bool steady_state = false;	//set to true once steady state is reached
	double last_mass = 0;	//mass at the prior time step
	double last_mom = 0;	//momentum at the prior time step
	double last_en = 0;		//energy at the prior time step

	std::chrono::time_point<std::chrono::high_resolution_clock> time_start;	//time at simulation start
};

#endif
