#ifndef SOURCE_H_
#define SOURCE_H_

#include "World.h"
#include "Species.h"

/*base class*/
class Source {
public:
	virtual void sample() = 0;
	virtual ~Source() {};	//to support destruction through base class
};


//simple monoenergetic source
class ColdBeamSource: public Source {
public:
	ColdBeamSource(Species &species, World &world, double v_drift, double den) :
		sp{species}, world{world}, v_drift{v_drift}, den{den} {}

	//generates particles
	void sample();

protected:
	Species &sp;	//reference to the injected species
	World &world;		//reference to world
	double v_drift;		//mean drift velocity
	double den;			//injection density
};

//simple monoenergetic source
class WarmBeamSource: public Source {
public:
	WarmBeamSource(Species &species, World &world, double v_drift, double den, double T) :
		sp{species}, world{world}, v_drift{v_drift}, den{den}, T{T} {}

	//generates particles
	void sample();

protected:
	Species &sp;	//reference to the injected species
	World &world;		//reference to world
	double v_drift;		//mean drift velocity
	double den;			//injection density
	double T;			//temperature
};


#endif /* SOURCE_H_ */
