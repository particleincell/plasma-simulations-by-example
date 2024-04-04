#ifndef SOURCE_H_
#define SOURCE_H_

#include "World.h"
#include "Species.h"

//simple monoenergetic source
class ColdBeamSource {
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


#endif /* SOURCE_H_ */
