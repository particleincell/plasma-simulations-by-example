/*definitions for species functions*/
#include <math.h>
#include "Species.h"
#include "Field.h"
#include <iostream>

/*loads randomly distributed particles in a x1-x2 box representing num_den number density*/
void Species::loadParticlesBox(double3 x1, double3 x2, double num_den, int num_mp)
{
	double box_vol = (x2[0]-x1[0])*(x2[1]-x1[1])*(x2[2]-x1[2]);		//box volume
	double num_real = num_den * box_vol;		//number of real particles
	double mpw = num_real/num_mp;			//macroparticle weight

	/*preallocate memory for particles*/
	particles.reserve(num_mp);

	/*load particles on an equally spaced grid*/
	for (int p=0;p<num_mp;p++)
	{
		//sample random position
		double3 pos;
		pos[0] = x1[0] + rnd()*(x2[0]-x1[0]);
		pos[1] = x1[1] + rnd()*(x2[1]-x1[1]);
		pos[2] = x1[2] + rnd()*(x2[2]-x1[2]);

		//set initial velocity
		double3 vel {0,0,0};	//stationary particle

		addParticle(pos,vel,mpw);	//add a new particle to the array
	}
}

/*quiet start load of num_sim[0]*num_sim[1]*num_sim[2] particles in a x1-x2 box
representing num_den number density*/
void Species::loadParticlesBoxQS(double3 x1, double3 x2, double num_den, int3 num_mp)
{
	double box_vol = (x2[0]-x1[0])*(x2[1]-x1[1])*(x2[2]-x1[2]);		//box volume
	int num_mp_tot = (num_mp[0]-1)*(num_mp[1]-1)*(num_mp[2]-1);	//total number of simulation particles
	double num_real = num_den * box_vol;		//number of real particles
	double mpw = num_real/num_mp_tot;			//macroparticle weight

	/*compute particle grid spacing*/
	double di = (x2[0]-x1[0])/(num_mp[0]-1);
	double dj = (x2[1]-x1[1])/(num_mp[1]-1);
	double dk = (x2[2]-x1[2])/(num_mp[2]-1);

	/*preallocate memory for particles*/
	particles.reserve(num_mp_tot);

	/*load particles on a equally spaced grid*/
	for (int i=0;i<num_mp[0];i++)
		for (int j=0;j<num_mp[1];j++)
			for (int k=0;k<num_mp[2];k++)
			{
				double pos[3];
				pos[0] = x1[0] + i*di;
				pos[1] = x1[1] + j*dj;
				pos[2] = x1[2] + k*dk;

				//shift particles on max faces back to the domain
				if (pos[0]==x2[0]) pos[0]-=1e-4*di;
				if (pos[1]==x2[1]) pos[1]-=1e-4*dj;
				if (pos[2]==x2[2]) pos[2]-=1e-4*dk;

				double w = 1;	//relative weight
				if (i==0 || i==num_mp[0]-1) w*=0.5;
				if (j==0 || j==num_mp[1]-1) w*=0.5;
				if (k==0 || k==num_mp[2]-1) w*=0.5;
				
				/*add rewind*/
				double vel[3] = {0,0,0};	//particle is stationary

				addParticle(pos,vel,mpw*w);	//add a new particle to the array
			}
}

/*adds a new particle, rewinding velocity by half dt*/
void Species::addParticle(double3 pos, double3 vel, double mpw)
{
    //add to list
    particles.emplace_back(pos,vel,mpw);
}
	
/*compute number density*/
void Species::computeNumberDensity()
{
	den.clear();		//set density to zero
	for (Particle &part:particles)
	{
		double3 lc = world.XtoL(part.pos);
		den.scatter(lc, part.mpw);
	}
		
	//divide by node volume
	den /= world.node_vol;
}



