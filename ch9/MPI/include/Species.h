/*Defines flying material data*/

#ifndef _SPECIES_H
#define _SPECIES_H

#include <list>
#include "World.h"
#include "Field.h"
#include <stdexcept>

using namespace std;

/** Data structures for particle storage **/
struct Particle
{
	double pos[3];			/*position*/
	double vel[3];			/*velocity*/
	int id  = -1;					/*particle id for output*/
	bool alive = false;				/*flag to mark dead particles*/

	Particle() {/*do nothing*/}	/*this gets called on memory allocation*/

	/*moved from constructor to "init" function called on add*/
	void init(double x[3], double v[3], int id)
		{for (int i=0;i<3;i++) {pos[i]=x[i];vel[i]=v[i];} this->id=id;alive=true;}
};

/*species container*/
class Species 
{
public:
	Species(World &world,string name, double mass, double charge, double spwt, int to_allocate) :
		max_allocated(to_allocate),vel(world),vel_ave(world),world(world) {

		this->name = name;
		this->mass=mass;
		this->charge=charge;
		this->spwt=spwt;
		den = new Field(world); 
		den_ave = new Field(world);

		particles = new Particle[max_allocated];
		if (!particles) throw runtime_error("Insufficient memory");
	}

	~Species() {
		delete den;delete den_ave;
		delete[] particles;
	}

	/*adds new particle*/
	void addParticle(double x[3], double v[3], int id=-1) {
		if (np<max_allocated)
			particles[np++].init(x,v,id>=0?id:part_id++);
	}

	/*moves all particles using electric field ef[]*/
	void move();

	/*moves particles between [start,end)*/
	static void moveKernel (Species *sp, int start, int end);

	/*sets specific weight*/
	void setSpwt(double spwt) {this->spwt = spwt;}	
	
	/*returns species total momentum*/
	double getMomentum();

	/*computes number density and velocity*/
	void computeGasProperties();

/*variables*/
	double mass;			/*particle mass in kg*/
	double charge;			/*particle charge in Coulomb*/
	double spwt;			/*species specific weight*/
	string name;
	
	/*particles now stored in array*/
	Particle *particles;
	int np = 0;
	const int max_allocated = 0;

	Field *den;
	Field *den_ave;		/*average density*/
	Field3 vel;
	Field3 vel_ave;

protected:
	World &world;
	int part_id = 0;

	/*MPI particle transfer*/
	list<Particle*> transfer[6];
	void transferParticles();

};


#endif
