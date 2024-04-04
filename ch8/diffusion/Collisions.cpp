#include <vector>
#include <math.h>
#include <iostream>
#include "Collisions.h"

using namespace std;


/*contant cross-section for the demo*/
double evalSigma(double rel_g)
{
	return 1e-16;
}


/*MCC*/
void MCC_CEX::apply(double dt)
{
		/*set pointers to target data*/
		Field &target_den = target.den;
		Field &target_temp = target.T;
		Field3 &target_vel = target.vel;

		/*loop over all particles*/
        for (Particle &part:source.particles)
        {
			/*get target velocity*/
			double3 lc = world.XtoL(part.pos);
            			
			double3 vt = target_vel.gather(lc);

            /*get target density*/
            double nn = target_den.gather(lc);
           
			/*compute cross-section, function of relative velocity*/
            double3 v_rel = part.vel-vt;
			double v_rel_mag = mag(v_rel);
            
            /*evaluate cross-section */
            double sigma = 1e-16;
            
            /*compute probability*/
            double P = 1 - exp(-nn*sigma*v_rel_mag*dt);
            
            /*compare to a random number to see if collision happened*/
            if (P>=rnd())
            {
                /*sample a virtual target particle*/
                double T_target = target_temp.gather(lc);	//sample target temperature
                double3 v_target = target.sampleIsotropicVel(T_target);
                part.vel = v_target; 	//CEX collision
            }
        }
    }
    

