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
		den(world.nn),	den_old(world.nn), den_ave(world.nn), vel(world.nn), j(world.nn),	world(world) { 	}

	/*returns the number of simulation particles*/
	size_t getNp()	{return particles.size();}

	/*returns the number of real particles*/
	double getRealCount();

	/*returns the species momentum*/
	double3 getMomentum();

	/*returns the species kinetic energy*/
	double getKE();

	/*moves all particles using electric field ef[]*/
	void advance();

	/*compute number density*/
	void computeGasProperties();

	/*samples velocity moments*/
	void sampleMoments();

	/*adds a new particle*/
	void addParticle(double3 pos, double3 vel) {addParticle(pos,vel,world.getDt());}
	void addParticle(double3 pos, double3 vel, double dt) {addParticle(pos,vel,dt,mpw0);}
	void addParticle(double3 pos, double3 vel, double dt, double mpw);

	/*updates number density*/
	void updateAverages() {den_ave.updateAverage(den);}

	/*returns random thermal velocity*/
	double sampleVth(double T);

	const std::string name;			/*species name*/
	const double mass;			/*particle mass in kg*/
	const double charge;		/*particle charge in Coulomb*/
	const double mpw0;			/*default macroparticle weight*/
	
	std::vector<Particle> particles;	/*contiguous array for storing particles*/
	Field den;			//number density
	Field den_old;		//prior time step density
	Field den_ave;		//averaged number density
	Field3 vel;			//velocity
	Field3 j;			//current density

protected:
	World &world;

};


#endif
