/*definitions for species functions*/
#include <math.h>
#include <iostream>
#include "Species.h"
#include "Field.h"

/*updates velocities and positions of all particles of this species*/
void Species::advance(Species &neutrals, Species &spherium)
{
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

			/*did this particle leave the domain?*/
			if (!world.inBounds(part.pos))
			{
				part.mpw = 0;	//kill the particle
			}
			else if (world.inSphere(part.pos)) {
				double tp = world.lineSphereIntersect(pos_old,part.pos);
				double dt_rem = (1-tp)*part.dt;
				part.dt -= dt_rem;

				//move particle *almost* to the surface
				part.pos = 	pos_old +0.999*tp*(part.pos-pos_old);
				double v_mag1 = mag(part.vel);	//pre-impact speed
				if (charge==0)	/*neutrals*/
					part.vel = sampleReflectedVelocity(part.pos,v_mag1);
				else {	/*ions*/
					double mpw_ratio = this->mpw0/neutrals.mpw0;
					part.mpw = 0;	//kill source particle

					//inject neutrals
					int mp_create = (int)(mpw_ratio+rnd());
					for (int i=0;i<mp_create;i++) {
						double3 vel = sampleReflectedVelocity(part.pos,v_mag1);
						neutrals.addParticle(part.pos,vel);
					}

					//emit sputtered material using a simple yield model
					double sput_yield = (v_mag1>5000)?0.1:0;
					double sput_mpw_ratio  = sput_yield*this->mpw0/spherium.mpw0;
					int sput_mp_create = (int)(sput_mpw_ratio+rnd());
					for (int i=0;i<sput_mp_create;i++) {
						double3 vel = sampleReflectedVelocity(part.pos,v_mag1);
						spherium.addParticle(part.pos,vel);
					}
				}
				continue;
			}

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
	return 3/sqrt(2+2+2)*sqrt(v1*v1+v2*v2+v3*v3);	//magnitude
}

/*returns random isotropic velocity*/
double3 Species::sampleIsotropicVel(double T) {
	double theta = 2*Const::PI*rnd();
    double r = -1.0+2*rnd();    //pick a random direction for d[0]
    double a = sqrt(1-r*r); //scaling for unity magnitude

    double3 d;
    d[0] = r;
	d[1] = cos(theta)*a;
    d[2] = sin(theta)*a;

    double v_th = sampleVth(T);
    double3 vel = v_th*d;
    return vel;
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
				double count = n_sum(i,j,k);
				if (count<=0) {T[i][j][k] = 0; continue;}

				double u_ave = vel(i,j,k)[0];
				double v_ave = vel(i,j,k)[1];
				double w_ave = vel(i,j,k)[2];
				double u2_ave = nuu_sum(i,j,k)/count;
				double v2_ave = nvv_sum(i,j,k)/count;
				double w2_ave = nww_sum(i,j,k)/count;

				double uu = u2_ave - u_ave*u_ave;
				double vv = v2_ave - v_ave*v_ave;
				double ww = w2_ave - w_ave*w_ave;
				T[i][j][k] = mass/(2*Const::K)*(uu+vv+ww);
			}
}

/*computes number of macroparticles per cell*/
void Species::computeMPC() {
	mpc.clear();
	for (Particle &part:particles) {
		int3 ijk = world.XtoIJK(part.pos);
		int i = ijk[0], j = ijk[1], k=ijk[2];
		mpc[i][j][k] += 1;
	}
}

/*clears sampled moment data*/
void Species::clearSamples() {
	n_sum = 0; nv_sum = 0; nuu_sum = 0; nvv_sum = 0; nww_sum=0;
}

