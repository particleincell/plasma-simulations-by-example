/*definitions for species functions*/
#include <math.h>
#include <iostream>
#include "Species.h"
#include "Field.h"

/*updates velocities and positions of all particles of this species*/
void Species::advance()
{
	double3 x0 = world.getX0();
	double3 xm = world.getXm();
	double3 dh = world.getDh();
	double Ly = dh[1]*(world.nj-1);

	/*loop over all particles*/
	for (Particle &part: particles)
	{
		/*increment particle's dt by world dt*/
		part.dt += world.getDt();

		/*get logical coordinate of particle's position*/
		double3 lc = world.XtoL(part.pos);
		
		/*electric field at particle position*/
		double3 ef_part = world.ef.gather(lc);
			
		/*update velocity from F=qE*/
		part.vel += ef_part*(part.dt*charge/mass);

		/*update position from v=dx/dt, take into account particle bounces*/
		int n_bounces = 0;

		/*keep iterate while time remains and the particle is alive*/
		while (part.dt>0 && part.mpw>0) {
			double3 pos_old = part.pos;
			part.pos += part.vel*part.dt;

			/*did the particle leave through the periodic boundary?*/
			if (part.pos[1]<x0[1]) 	part.pos[1] += Ly;
			else if (part.pos[1]>=xm[1]) part.pos[1] -= Ly;
			
			/*did this particle leave the domain or hit the sphere?*/
			if (!world.inBounds(part.pos) || 
				 world.inSphere(part.pos))	part.mpw = 0;	//kill the particle
			
			//this particle finished the whole step
			part.dt = 0;

			//kill stuck particles
			if (++n_bounces>20) {std::cerr<<"Stuck particle!"<<std::endl;part.mpw = 0;}
		}
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

/*returns random post-impact velocity*/
double3 Species::sampleReflectedVelocity(double3 pos, double v_mag1)
{
	double v_th = sampleVth(1000); //assume T_sphere = 1000K
	const double a_th = 1;		//thermal accommodation coeff
	double v_mag2 = v_mag1 + a_th*(v_th-v_mag1);
	return v_mag2*world.sphereDiffuseVector(pos); //set new velocity
}

/*adds a new particle, rewinding velocity by half dt*/
void Species::addParticle(double3 pos, double3 vel, double dt, double mpw)
{
	//don't do anything (return) if pos outside domain bounds [x0,xd)
	if (!world.inBounds(pos)) return;

	//get particle logical coordinate
	double3 lc = world.XtoL(pos);
	
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

/*returns random thermal velocity*/
double Species::sampleVth(double T)
{
	//thermal velocity
	double v_th = sqrt(2*Const::K*T/mass);
	//get three random velocity components
	double v1 = v_th*(rnd()+rnd()+rnd()-1.5);
	double v2 = v_th*(rnd()+rnd()+rnd()-1.5);
	double v3 = v_th*(rnd()+rnd()+rnd()-1.5);
	return sqrt(v1*v1+v2*v2+v3*v3);	//magnitude
}

/*compute number density*/
void Species::computeNumberDensity()
{
	den.clear();
	for (Particle &part:particles)
	{
		double3 lc = world.XtoL(part.pos);
		den.scatter(lc, part.mpw);
	}

	/*add up contributions along the periodic boundary*/
	for (int i=0;i<world.ni;i++)
		for (int k=0;k<world.nk;k++) {
			den[i][0][k] += den[i][world.nj-1][k];
			den[i][world.nj-1][k] = den[i][0][k];
	}


	//divide by node volume
	den /= world.node_vol;
}

/*samples velocity moments*/
void Species::sampleMoments() {
	for (Particle &part:particles)
	{
		double3 lc = world.XtoL(part.pos);
		n_sum.scatter(lc, part.mpw);
		nv_sum.scatter(lc,part.mpw*part.vel);
		nuu_sum.scatter(lc,part.mpw*part.vel[0]*part.vel[0]);
		nvv_sum.scatter(lc,part.mpw*part.vel[1]*part.vel[1]);
		nww_sum.scatter(lc,part.mpw*part.vel[2]*part.vel[2]);
	}
}

/*uses sampled data to compute velocity and temperature*/
void Species::computeGasProperties() {
	vel = nv_sum/n_sum;	//stream velocity

	for (int i=0;i<world.ni;i++)
		for (int j=0;j<world.nj;j++)
			for (int k=0;k<world.nk;k++) {
				int count = n_sum(i,j,k);
				double u_ave = vel(i,j,k)[0];
				double v_ave = vel(i,j,k)[1];
				double w_ave = vel(i,j,k)[2];
				double u2_ave = (count>0)?(nuu_sum(i,j,k)/count):0;
				double v2_ave = (count>0)?(nvv_sum(i,j,k)/count):0;
				double w2_ave = (count>0)?(nww_sum(i,j,k)/count):0;

				double uu = u2_ave - u_ave*u_ave;
				double vv = v2_ave - v_ave*v_ave;
				double ww = w2_ave - w_ave*w_ave;
				if (!isfinite(uu+vv+ww))
				{
					std::cout<<"infinite T"<<std::endl;
				}
				T[i][j][k] = mass/(2*Const::K)*(uu+vv+ww);
			}
}

/*clears sampled moment data*/
void Species::clearSamples() {
	n_sum = 0; nv_sum = 0; nuu_sum = 0; nvv_sum = 0; nww_sum=0;
}

