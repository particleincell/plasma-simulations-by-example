/*
Definition for a particle source
 */

#include "Species.h"
#include "World.h"
#include <math.h>

class Source 
{
public:

	World &world;
    Species &species;    /*injection species*/
    double den0;        /*injection density*/
    double v_drift;     /*source velocity*/
    double v_th;		/*thermal velocity*/
    
    Source(World &world, Species &species, double den0, double v_drift, double T):
		world(world), species(species)
    {
        this->den0 = den0;
        this->v_drift = v_drift;
        this->v_th = sqrt(2*Kb*T/species.mass);
		num_rem = 0;
    }
    
    
    /*samples ions at the k=0 edge such that local plasma density equals PLASMA_DEN*/
	double num_rem;		/*remainder number of ions*/	
	void sample();
};
