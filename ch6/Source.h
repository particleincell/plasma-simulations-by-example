#ifndef SOURCE_H_
#define SOURCE_H_

#include "World.h"
#include "Species.h"

//simple monoenergetic source
class ColdBeamSource {
public:
	ColdBeamSource(Species &species, World &world, double v_drift, double den, Group &src_group);

	//generates particles
	void sample();

protected:
	Species &sp;	//reference to the injected species
	World &world;		//reference to world
	double v_drift;		//mean drift velocity
	double den;			//injection density
	double A_tot;
	Group &src_group;	//injection surface triangles
};


#endif /* SOURCE_H_ */
