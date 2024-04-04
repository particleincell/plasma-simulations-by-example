#include "Source.h"


//samples monoenergetic particles according to a prescribed density
void ColdBeamSource::sample()
{
	double3 dh = world.getDh();
	double3 x0 = world.getX0();

	//area of the XY plane, A=Lx*Ly
	double Lx = dh[0]*(world.ni-1);
	double Ly = dh[1]*(world.nj-1);
	double A = Lx*Ly;

	//compute number of real particles to generate: (#/s) = n*v*A; # = (#/s)*dt
	double num_real = den*v_drift*A*world.getDt();

	//number of simulation particles
	int num_sim = (int)(num_real/sp.mpw0+rnd());

	//inject particles
	for (int i=0;i<num_sim;i++)
	{
		double3 pos {x0[0]+rnd()*Lx, x0[1]+rnd()*Ly, x0[2]};
		double3 vel {0,0,v_drift};
		sp.addParticle(pos,vel);
	}
}

//samples particles with finite thermal and drift velocity
void WarmBeamSource::sample()
{
	double3 dh = world.getDh();
	double3 x0 = world.getX0();

	//area of the XY plane, A=Lx*Ly
	double Lx = dh[0]*(world.ni-1);
	double Ly = dh[1]*(world.nj-1);
	double A = Lx*Ly;

	//compute number of real particles to generate: (#/s) = n*v*A; # = (#/s)*dt
	double num_real = den*v_drift*A*world.getDt();

	//number of simulation particles
	int num_sim = (int)(num_real/sp.mpw0+rnd());

	//inject particles
	for (int i=0;i<num_sim;i++)
	{
		double3 pos {x0[0]+rnd()*Lx, x0[1]+rnd()*Ly, x0[2]};
		double3 vel = sp.sampleIsotropicVel(T);
        vel[2] += v_drift; //add drift component

        sp.addParticle(pos,vel);
	}
}
