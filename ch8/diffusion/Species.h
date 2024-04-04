/*Defines flying material data*/

#ifndef _SPECIES_H
#define _SPECIES_H

#include <vector>
#include <sstream>
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
class Species {
public:

	Species(std::string name, double mass, double charge, World &world) :
		name(name), mass(mass), charge(charge),
		den(world.nn), T(world.nn),	vel(world.nn),
		den_ave(world.nn),
		world(world) { 	}
	~Species() {};


	virtual void advance() = 0; //integrates species state by dt

	double3 sampleIsotropicVel(double T);	//samples random isotropic velocity
	double sampleVth(double T);		//returns random thermal velocity

	/*updates number density*/
	void updateAverages() {den_ave.updateAverage(den);}

	/*hooks to be overridden as needed*/
	virtual void computeGasProperties() {} // uses sampled data to compute velocity and temperature
	virtual void clearSamples() {} //clears averaged data

	virtual double getMass() {return 0;}		//total mass
	virtual double3 getMomentum() {return {0,0,0};}	//total momentum
	virtual double getKE() {return 0;}			//total kinetic energy


	virtual std::string printSelf() = 0;

	const std::string name;			/*species name*/
	const double mass;			/*particle mass in kg*/
	const double charge;		/*particle charge in Coulomb*/

	Field den;			/*number density*/
	Field T;			/*temperature*/
	Field3 vel;			/*stream velocity*/
	Field den_ave;		/*averaged number density*/

protected:
	World &world;

};

class FluidSpecies:public Species {
public:
	FluidSpecies(std::string name, double mass, double charge, double den0, double D,World &world):
		Species{name, mass, charge, world}, den0{den0}, D{D}, den_new(world.nn) {}

	void advance();
	std::string printSelf();

protected:
	double den0;	//sphere density
	double D;		//diffusion coefficient
	Field den_new;	//temporary solution buffer
};

class KineticSpecies:public Species
{
public:
	KineticSpecies(std::string name, double mass, double charge, double mpw0, World &world) :
		Species(name, mass, charge, world),
		 mpw0(mpw0), mpc(world.ni-1,world.nj-1,world.nk-1),
		n_sum(world.nn),nv_sum(world.nn),
		nuu_sum(world.nn), nvv_sum(world.nn), nww_sum(world.nn)		{ 	}

	std::string printSelf();


	/*returns the number of simulation particles*/
	size_t getNp()	{return particles.size();}

	/*returns the number of real particles*/
	double getMass();

	/*returns the species momentum*/
	double3 getMomentum();

	/*returns the species kinetic energy*/
	double getKE();

	/*moves all particles using electric field ef[]*/
	void advance();

	/*compute number density*/
	void computeNumberDensity();

	/*samples velocity moments*/
	void sampleMoments();

	/*uses sampled data to compute velocity and temperature*/
	void computeGasProperties();

	/*computes number of macroparticles per cell*/
	void computeMPC();

	double3 sampleReflectedVelocity(double3 pos, double v_mag1);  //returns random post-impact velocity*/


	/*adds a new particle*/
	void addParticle(double3 pos, double3 vel) {addParticle(pos,vel,world.getDt());}
	void addParticle(double3 pos, double3 vel, double dt) {addParticle(pos,vel,dt,mpw0);}
	void addParticle(double3 pos, double3 vel, double dt, double mpw);

	/*clears sampled moment data*/
	void clearSamples();


	const double mpw0;			/*default macroparticle weight*/
	
	std::vector<Particle> particles;	/*contiguous array for storing particles*/
	Field mpc;			/*macroparticles per cell*/

protected:

	Field n_sum;
	Field3 nv_sum;
	Field nuu_sum,nvv_sum,nww_sum;

};



#endif
