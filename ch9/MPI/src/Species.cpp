/*definitions for species functions*/
#include "Species.h"
#include "Field.h"
#include "PotentialSolver.h"	/*for vec*/
#include <thread>
#include <mpi.h>

namespace vec {
	void unit(double r[3], double v1[3])
	{
		double v_mag = sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
		
		for (int i=0;i<3;i++) r[i]=v1[i]/v_mag;	
	}
    
	void add(double r[3], double v1[3], double v2[3])
	{    
		for (int i=0;i<3;i++) r[i]=v1[i]+v2[i];
	}
    
	void sub(double r[3], double v1[3], double v2[3])
	{    
		for (int i=0;i<3;i++) r[i]=v1[i]-v2[i];
	}
    
	double dot(double v1[3], double v2[3])
	{
		double r=0;
		for (int i=0;i<3;i++) r+=v1[i]*v2[i];
		return r;
	}
    
	void cross(double r[3], double v1[3], double v2[3])
	{
		r[0] = v1[1]*v2[2]-v1[2]*v2[1];
		r[1] = -v1[0]*v2[2]+v1[2]*v2[0];
		r[2] = v1[0]*v2[1]-v1[1]*v2[0];
	}
}   

void Species::moveKernel(Species *sp, int start, int end)
{
	/*pointer to ef*/
	World &world = sp->world;
	Field3 *ef = world.ef;
	double qm = sp->charge/sp->mass;

	/*continue while particles remain*/
	for (int p=start;p<end;p++)
	{
		Particle &part = sp->particles[p];
		
		/*get logical coordinate of particle's position*/
		double lc[3];
		world.XtoL(lc, part.pos);
		
		/*electric field at particle position*/
		double ef_part[3];
		ef->gather(ef_part,lc);
			
		/*update velocity from F=qE*/
		for (int i=0;i<3;i++)
			part.vel[i] += qm*ef_part[i]*world.dt;
			
		double x_old[3];		/*position before the push*/
		for (int i=0;i<3;i++) x_old[i] = part.pos[i];

		/*update position from v=dx/dt*/
		for (int i=0;i<3;i++)
			part.pos[i] += part.vel[i]*world.dt;
			

		/*first check for sphere hit, it may result in new particle position if reflecting back*/
		if (world.inSphere(part.pos))	/*did this particle get inside the sphere?*/
		{
			if (sp->charge!=0) part.alive = false;   /*remove ions, these should be turned into neutrals!*/
			else
			{
			   /*diffusely reflect neutrals*/

				/*approximate intersection point, move particle there*/
				for (int i=0;i<3;i++)
					part.pos[i] = 0.5*(x_old[i]+part.pos[i]);

				/*compute normal and two tangents*/
				double r[3],norm[3];
				vec::sub(r,part.pos,world.sphere_x0);
				vec::unit(norm,r);

				/*two perpendicular vectors,
				http://math.stackexchange.com/questions/137362/how-to-find-perpendicular-vector-to-another-vector*/
				double perp1[3] = {norm[1],-norm[0],0};
				double perp2[3] = {-norm[2],0,norm[0]};
				double tang1[3];

				vec::add(tang1,perp1,perp2);
				vec::unit(tang1,tang1);

				/*final tang given by cross product*/
				double tang2[3];
				vec::cross(tang2,norm,tang1);
				vec::unit(tang2,tang2);

				/*sanity check*/
				/*
				double dot1 = vec_dot(norm,tang1);
				double dot2 = vec_dot(norm,tang2);
				double dot3 = vec_dot(tang2,tang1);
				*/

				/*sample velocity from cosine distribution*/
				double sin_theta = rnd();
				double cos_theta = sqrt(1-sin_theta*sin_theta);

				//random in plane angle
				double psi = rnd()*2*PI;

				//three vector components
				double a = sin_theta*cos(psi);
				double b = sin_theta*sin(psi);
				double c = cos_theta;

				//multiply by corresponding directions
				double v1[3];
				double v2[3];
				double v3[3];

				for (int i=0;i<3;i++)
				{
					v1[i] = tang1[i]*a;
					v2[i] = tang2[i]*b;
					v3[i] = norm[i]*c;
				}

				 /*get initial velocity magnitude*/
				double v_mag = sqrt(part.vel[0]*part.vel[0]+
						   part.vel[1]*part.vel[1]+part.vel[2]*part.vel[2]);

				/*add velocity components*/
				for (int i=0;i<3;i++) part.vel[i] = v_mag*(v1[i]+v2[i]+v3[i]);
			}
		} /*if in sphere*/

		/*if still alive (not absorbed by the sphere*/
		if (part.alive)
		{
			/*compute new logical coordinate (slow!) for boundary checking
			  could instead compare against world.x0 / world.xd	*/
			world.XtoL(lc, part.pos);

			/*check if particle crossed processor boundary*/
			if (world.mpi_i>0 && lc[0]<0)
				sp->transfer[0].emplace_back(&part);
			else if (world.mpi_i<world.mpi_size_i-1 && lc[0]>=world.ni-1)
				sp->transfer[1].emplace_back(&part);
			else if (world.mpi_j>0 && lc[1]<0)
				sp->transfer[2].emplace_back(&part);
			else if (world.mpi_j<world.mpi_size_j-1 && lc[1]>=world.nj-1)
				sp->transfer[3].emplace_back(&part);
			else if (world.mpi_k>0 && lc[2]<0)
				sp->transfer[4].emplace_back(&part);
			else if (world.mpi_k<world.mpi_size_k-1 && lc[2]>=world.nk-1)
				sp->transfer[5].emplace_back(&part);

			/*kill particle if leaving the domain (this includes crossing MPI boundary)?*/
			if(lc[0]<0 || lc[0]>=world.ni-1 ||
			   lc[1]<0 || lc[1]>=world.nj-1 ||
			   lc[2]<0 || lc[2]>=world.nk-1)
			   part.alive = false;
		}

	} //for p

}

#define min(a,b) (a<b?a:b)
/*updates velocities and positions of all particles of this species*/
void Species::move()
{
	if (world.num_threads>1)
	{
		/*multithreaded launch*/
		/*thread chunk size*/
		int chunk = np/world.num_threads+1;
		thread *threads = new thread[world.num_threads];

		for (int i=0;i<world.num_threads;i++)
			threads[i] = thread(moveKernel,this, i*chunk,min((i+1)*chunk,np));

		for (int i=0;i<world.num_threads;i++)
			threads[i].join();
		delete[] threads;
	}
	else
		moveKernel(this,0,np);	/*single threaded*/

	/*transfer particles between MPI processes*/
	transferParticles();

	/*perform particle removal sweep*/
	for (int p=0;p<np;p++)
	{
		if (!particles[p].alive)
			particles[p--] = particles[--np];
	}

	/*compute density*/
	computeGasProperties();
}

void Species::transferParticles()
{
	/**** PARTICLE TRANSFER ****/
	const int TAG_COUNT = 20;		/*some tag*/
	const int TAG_DATA = 22;


	for (int face=0;face<6;face++)
	{
		int target;

		target = world.getNeighborRank(face);

		/*transfer if we have a neighbor*/
		if (target>=0)
		{
			/*start by transferring particle counts*/
			int send_size = transfer[face].size();
			int recv_size = 0;

			MPI_Sendrecv(&send_size, 1, MPI_INT, target, TAG_COUNT,
						 &recv_size, 1, MPI_INT, target, TAG_COUNT,
						 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			/*allocate memory buffers for passing particles*/
			double *send_dbuffer=nullptr, *recv_dbuffer=nullptr;
			int *send_ibuffer=nullptr, *recv_ibuffer=nullptr;

			/*since we can't allocate zero memory*/
			if (send_size>0)
			{
				send_dbuffer = new double[send_size*6];	/*6 doubles per particle*/
				send_ibuffer = new int[send_size*1];	/*1 int per particle*/
			}
			if (recv_size>0)
			{
				recv_dbuffer = new double[recv_size*6];
				recv_ibuffer = new int[recv_size*1];
			}

			/*pack particles to a buffer*/
			auto it = transfer[face].begin();
			for (int i = 0;i<send_size;i++,it++)
			{
				Particle &part = **it; // *it is pointer so **it is the reference
				send_dbuffer[6*i+0] = part.pos[0];
				send_dbuffer[6*i+1] = part.pos[1];
				send_dbuffer[6*i+2] = part.pos[2];
				send_dbuffer[6*i+3] = part.vel[0];
				send_dbuffer[6*i+4] = part.vel[1];
				send_dbuffer[6*i+5] = part.vel[2];
				send_ibuffer[i] = part.id;
			}
			/*transfer particles*/
			MPI_Sendrecv(send_dbuffer, send_size*6, MPI_DOUBLE, target, TAG_DATA,
						 recv_dbuffer, recv_size*6, MPI_DOUBLE, target, TAG_DATA,
						 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Sendrecv(send_ibuffer, send_size*1, MPI_INT, target, TAG_DATA,
						 recv_ibuffer, recv_size*1, MPI_INT, target, TAG_DATA,
						 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			/*append particles from received buffer to our particle list*/
			for (int i = 0;i<recv_size;i++)
			{
				double pos[3];
				double vel[3];

				pos[0] = recv_dbuffer[6*i+0];
				pos[1] = recv_dbuffer[6*i+1];
				pos[2] = recv_dbuffer[6*i+2];
				vel[0] = recv_dbuffer[6*i+3];
				vel[1] = recv_dbuffer[6*i+4];
				vel[2] = recv_dbuffer[6*i+5];
				int id = recv_ibuffer[i];

				/*only add particle if it belongs to us*/
				if (pos[0]>=world.x0[0] && pos[0]<world.xd[0] &&
					pos[1]>=world.x0[1] && pos[1]<world.xd[1] &&
					pos[2]>=world.x0[2] && pos[2]<world.xd[2])
						addParticle(pos,vel,id);
				else
				{
					/*discard particle for simplicity, not right!!*/
				}
			}

			/*free memory*/
			if (send_size>0)
			{
				delete[] send_dbuffer;
				delete[] send_ibuffer;
			}
			if (recv_size>0)
			{
				delete[] recv_dbuffer;
				delete[] recv_ibuffer;
			}
		} /*transfer*/

		/*empty transfer list*/
		transfer[face].clear();

	} /*for face*/
}

/*compute number density*/
void Species::computeGasProperties()
{
	den->clear();
	for (int p=0;p<np;p++)
	{
		Particle &part = particles[p];
		double lc[3];
		world.XtoL(lc, part.pos);
		den->scatter(lc, spwt);
		vel.scatter(lc,part.vel,spwt);
	}

	/*update data across MPI domain*/
	den->updateBoundaries();
	for (int i=0;i<3;i++) vel.f[i]->updateBoundaries();

	/*right now, density contains number of particles 
	 *per node, use this to get node-average velocity*/
    for (int i=0;i<3;i++)
        vel.f[i]->divideByField(den);

	den->divideByVolume();
}

/*computes total momentum*/
double Species::getMomentum()
{
	    double momentum  = 0;
        for (int p=0;p<np;p++)
        {
        	Particle &part = particles[p];
            double v_mag = sqrt(part.vel[0]*part.vel[0]+
                        part.vel[1]*part.vel[1]+part.vel[2]*part.vel[2]);
            momentum += v_mag*mass*spwt;            
        }
        return momentum;
}
