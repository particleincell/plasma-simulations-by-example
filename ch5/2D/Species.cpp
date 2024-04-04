/*definitions for species functions*/
#include <math.h>
#include <iostream>
#include "Species.h"
#include "Field.h"

/*updates velocities and positions of all particles of this species*/
void Species::advance()
{
	/*loop over all particles*/
	for (Particle &part: particles)
	{
		/*increment particle's dt by world dt*/
		part.dt += world.getDt();

		/*get logical coordinate of particle's position*/
		double2 lc = world.XtoL(part.pos);
		
		/*electric field at particle position*/
		double3 ef_part = world.ef.gather(lc);
			
		/*update velocity from F=qE*/
		part.vel += ef_part*(part.dt*charge/mass);

		/*update position from v=dx/dt, take into account particle bounces*/
		double3 pos_old = part.pos;
		part.pos += part.vel*part.dt;

		/*did this particle leave the domain?*/
		if (!world.inBounds(part.pos))
		{	
			part.mpw = 0;	//kill the particle
		}
		//check for particle hitting the object		
		else if (world.inCircle(part.pos)) {
			part.mpw = 0;	//kill source particle
		}
		
		//this particle finished the whole step
		part.dt = 0;
	}

	/*perform a particle removal step, dead particles are replaced by the entry at the end*/
	size_t np = particles.size();
	for (size_t p=0;p<np;p++)
	{
		if (particles[p].mpw>0) continue;	//ignore live particles
		particles[p] = particles[np-1]; //fill the hole
		np--;	//reduce count of valid elements
		p--;	//decrement p so this position gets checked again
	}

	//now delete particles[np:end]
	particles.erase(particles.begin()+np,particles.end());
}

/*adds a new particle, rewinding velocity by half dt*/
void Species::addParticle(double3 pos, double3 vel, double dt, double mpw)
{
	//don't do anything (return) if pos outside domain bounds [x0,xd)
	if (!world.inBounds(pos)) return;

	//get particle logical coordinate
	double2 lc = world.XtoL(pos);
	
	//evaluate electric field at particle position
    double3 ef_part = world.ef.gather(lc);

	//rewind velocity back by 0.5*dt*ef
    vel -=  charge/mass*ef_part*(0.5*world.getDt());

    //add to list
    particles.emplace_back(pos,vel,dt,mpw);
}

/*returns the number of real particles*/
double Species::getRealCount() {
	double mpw_sum = 0;
	for (Particle &part:particles)
		mpw_sum+=part.mpw;
	return mpw_sum;
}

/* returns the species momentum*/
double3 Species::getMomentum() {
	double3 mom;
	for (Particle &part:particles)
		mom+=part.mpw*part.vel;
	return mass*mom;
}

/* returns the species kinetic energy*/
double Species::getKE() {
	double ke = 0;
	for (Particle &part:particles)
	{
		double v2 = part.vel[0]*part.vel[0] + part.vel[1]*part.vel[1] + part.vel[2]*part.vel[2];
		ke += part.mpw*v2;
	}
	return 0.5*mass*ke;
}

/*compute number density*/
void Species::computeNumberDensity()
{
	den.clear();
	for (Particle &part:particles)
	{
		double2 lc = world.XtoL(part.pos);
		den.scatter(lc, part.mpw);
	}

	//divide by node volume
	den /= world.node_vol;
}


