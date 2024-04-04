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
	double mpw;				/*macroparticle weight*/
	LCord lc;				/*logical coordinate*/
	
	Particle(double3 x, double3 v, double mpw, LCord lc):
		pos{x}, vel{v}, mpw{mpw}, lc{lc} { }
};

class FESolver;

/*species container*/
class Species 
{
public:
	Species(std::string name, double mass, double charge, double mpw0, World &world) :
		name(name), mass(mass), charge(charge), mpw0(mpw0),
		den(world.vm), den_ave{world.vm}, world(world) { 	}

	/*returns the number of simulation particles*/
	size_t getNp()	{return particles.size();}

	/*returns the number of real particles*/
	double getRealCount();

	/*returns the species momentum*/
	double3 getMomentum();

	/*returns the species kinetic energy*/
	double getKE();

	/*moves all particles using electric field ef[]*/
	void advance(FESolver &solver);

	/*compute number density*/
	void computeNumberDensity();

	/*adds a new particle*/
	void addParticle(double3 pos, double3 vel, double mpwt);

	/*updates number density*/
	void updateAverages() {den_ave.updateAverage(den);}

	const std::string name;			/*species name*/
	const double mass;			/*particle mass in kg*/
	const double charge;		/*particle charge in Coulomb*/
	const double mpw0;			/*default macroparticle weight*/
	
	std::vector<Particle> particles;	/*contiguous array for storing particles*/
	Field den;			/*number density*/
	Field den_ave;		/*averaged number density*/

protected:
	World &world;

};


#endif
