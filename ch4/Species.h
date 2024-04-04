/*Defines flying material data*/

#ifndef _SPECIES_H
#define _SPECIES_H

#include <vector>
#include "Field.h"
#include "World.h"

/** Data structures for particle storage **/
struct Particle
{
	double3 pos;			/*position*/
	double3 vel;			/*velocity*/
	double dt;				/*time step to push the particle through*/
	double mpw;				/*macroparticle weight*/
	
	Particle(double3 x, double3 v, double dt, double mpw):
	pos{x}, vel{v}, dt{dt}, mpw{mpw} { }
};

/*species container*/
class Species 
{
public:
	Species(std::string name, double mass, double charge, double mpw0, World &world) :
		name(name), mass(mass), charge(charge), mpw0(mpw0),
		den(world.nn), T(world.nn),	vel(world.nn),
		den_ave(world.nn),mpc(world.ni-1,world.nj-1,world.nk-1),
		n_sum(world.nn),nv_sum(world.nn),
		nuu_sum(world.nn), nvv_sum(world.nn), nww_sum(world.nn),
		world(world) { 	}

	/*returns the number of simulation particles*/
	size_t getNp()	{return particles.size();}

	/*returns the number of real particles*/
	double getRealCount();

	/*returns the species momentum*/
	double3 getMomentum();

	/*returns the species kinetic energy*/
	double getKE();

	/*moves all particles using electric field ef[]*/
	void advance(Species &neutrals, Species &spherium);

	/*compute number density*/
	void computeNumberDensity();

	/*samples velocity moments*/
	void sampleMoments();

	/*uses sampled data to compute velocity and temperature*/
	void computeGasProperties();

	/*computes number of macroparticles per cell*/
	void computeMPC();

	/*clears sampled moment data*/
	void clearSamples();

	/*adds a new particle*/
	void addParticle(double3 pos, double3 vel) {addParticle(pos,vel,world.getDt());}
	void addParticle(double3 pos, double3 vel, double dt) {addParticle(pos,vel,dt,mpw0);}
	void addParticle(double3 pos, double3 vel, double dt, double mpw);

	/*updates number density*/
	void updateAverages() {den_ave.updateAverage(den);}

	/*returns random thermal velocity*/
	double sampleVth(double T);

	/*samples random isotropic velocity*/
	double3 sampleIsotropicVel(double T);

	const std::string name;			/*species name*/
	const double mass;			/*particle mass in kg*/
	const double charge;		/*particle charge in Coulomb*/
	const double mpw0;			/*default macroparticle weight*/
	
	std::vector<Particle> particles;	/*contiguous array for storing particles*/
	Field den;			/*number density*/
	Field T;			/*temperature*/
	Field3 vel;			/*stream velocity*/
	Field den_ave;		/*averaged number density*/
	Field mpc;			/*macroparticles per cell*/

protected:
	double3 sampleReflectedVelocity(double3 pos, double v_mag1);	/*returns random post-impact velocity*/

	Field n_sum;
	Field3 nv_sum;
	Field nuu_sum,nvv_sum,nww_sum;
	World &world;
};


#endif
