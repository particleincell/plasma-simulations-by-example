#include "Source.h"


//samples monoenergetic particles according to a prescribed density
void ColdBeamSource::sample()
{
	double3 dh = world.getDh();
	double3 x0 = world.getX0();

	//area of the sampling region plane, A=Lx*Ly
	double Ly = dh[1]*(0.4*world.nj-1); //sample small beam
	double A = Ly*1;

	//compute number of real particles to generate: (#/s) = n*v*A; # = (#/s)*dt
	double num_real = den*v_drift*A*world.getDt();

	//number of simulation particles
	int num_sim = (int)(num_real/sp.mpw0+rnd());

	//inject particles
	for (int i=0;i<num_sim;i++)
	{
		double3 pos {0.0*(world.ni-1)*dh[0], x0[1]+0.3*(world.nj-1)*dh[1]+rnd()*Ly, 0};
		double3 vel {v_drift,0,0};
		sp.addParticle(pos,vel);
	}
}

