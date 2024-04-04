#include "Source.h"


//samples monoenergetic particles according to a prescribed density
void ColdBeamSource::sample()
{
	double2 dh = world.getDh();
	double2 x0 = world.getX0();
	double2 xm = world.getXm();

	//area of a cirlce
	double Lr = xm[1]-x0[1];	//radial distance
	double A = Const::PI*(xm[1]*xm[1]-x0[1]*x0[1]);	//pi*(r2^2-r1^2)

	//compute number of real particles to generate: (#/s) = n*v*A; # = (#/s)*dt
	double num_real = den*v_drift*A*world.getDt();

	//number of simulation particles
	int num_sim = (int)(num_real/sp.mpw0+rnd());

	//inject particles
	for (int i=0;i<num_sim;i++)
	{
		double3 pos {x0[0], sqrt(Lr*Lr*rnd()+x0[1]*x0[1]), 0};
		double3 vel {v_drift,0,0};
		sp.addParticle(pos,vel);
	}
}

//P = 2*k*dr/(pi*r*r)
//R = r^2/(r_s^2)
//r = sqrt(R*r_s^2)
//r = sqrt(R)*r_s
