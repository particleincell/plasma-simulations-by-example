/*defines the simulation domain*/
#ifndef _WORLD_H
#define _WORLD_H

#include <math.h>
#include <iostream>
#include <sstream>

class Field;	/*tell the compiler that Field is a class that will be defined later*/
class Field3;

/*constants*/
const double EPS_0 = 8.85418782e-12;  	// F/m, vacuum permittivity
const double ME = 9.10938215e-31;		// kg, electron mass
const double QE = 1.602176565e-19;		// C, electron charge
const double AMU = 1.660538921e-27;		// kg, atomic mass unit
const double Kb = 1.3806488e-23;              // J/K, boltzmann constant
const double eVtoK = QE/Kb;                   // temperature conversion from eV to K
const double PI = acos(-1.0);

class World 
{
public:	
	/*FUNCTIONS*/

	/*constructor, allocats memory*/
	World(int ni, int nj, int nk);

	/*destructor, frees memory*/
	~World();

	/*functions to set mesh origin and spacing*/
	void setOrigin(double x0,double y0,double z0) {this->x0[0] = x0;this->x0[1]=y0;this->x0[2]=z0;}
	void setSpacing(double dx, double dy, double dz) {this->dh[0]=dx;this->dh[1]=dy;this->dh[2]=dz;}

	/*sets simulation time step*/
	void setDt(double dt) {this->dt=dt;}
	
	/*converts physical position to logical coordinate*/
	inline void XtoL(double lc[3], double pos[3])
	{
		for (int i=0;i<3;i++)
			lc[i] = (pos[i]-x0[i])/(dh[i]);
	}
	
	/*converts logical coordinate to physical position*/
	void pos(double x[3], double i, double j, double k)
	{
		double lc[3] = {i,j,k};
		pos(x, lc);
	}
	
	/*converts logical coordinate to physical position*/
	void pos(double x[3], double lc[3])
	{
		for (int i=0;i<3;i++)
			x[i] = x0[i]+dh[i]*lc[i];
	}

	/*marks k=0 plane as dirichlet*/
	void AddInlet();

	/*sugarcubes a sphere centered at (x0,y0,z0)*/
	void AddSphere(double x0, double y0, double z0, double radius, double phi_sphere);
	
	/*return true if point r is inside or on the sphere*/
	bool inSphere(double r[3]);

	/*sets number of threads*/
	void setNumberOfThreads(int num_threads) {this->num_threads = num_threads;}

	/*computes local domain from rank and domain decomposition*/
	void initMPIDomain(int mpi_rank, int mpi_size, int doms[3])
	{
		/*first save parameters*/
		this->mpi_rank = mpi_rank;
		this->mpi_size = mpi_size;
		this->mpi_size_i = doms[0];
		this->mpi_size_j = doms[1];
		this->mpi_size_k = doms[2];

		/*now set our position in the MPI world*/
		mpi_k = mpi_rank/(mpi_size_i*mpi_size_j);
		int mod = mpi_rank-mpi_k*(mpi_size_i*mpi_size_j);
		mpi_j = mod/mpi_size_i;
		mpi_i = mod%mpi_size_i;

		/*retain global nodes number*/
		num_nodes_global = ni*nj*nk;

		/*this is for clarity*/
		int ni_global = ni;
		int nj_global = nj;
		int nk_global = nk;

		/*number of nodes per process*/
		int proc_ni = ni/mpi_size_i;
		int proc_nj = nj/mpi_size_j;
		int proc_nk = nk/mpi_size_k;

		/*set the range we are simulating*/
		int my_i1 = proc_ni*mpi_i;
		int my_i2 = my_i1 + proc_ni+1;
		if (my_i2>ni_global) my_i2 = ni_global;

		int my_j1 = proc_nj*mpi_j;
		int my_j2 = my_j1 + proc_nj+1;
		if (my_j2>nj_global) my_j2 = nj_global;

		int my_k1 = proc_nk*mpi_k;
		int my_k2 = my_k1 + proc_nk+1;
		if (my_k2>nk_global) my_k2 = nk_global;

		/*set extents of our domain*/
		x0[0] = x0[0] + my_i1*dh[0];
		x0[1] = x0[1] + my_j1*dh[1];
		x0[2] = x0[2] + my_k1*dh[2];

		/*domain extent*/
		xd[0] = x0[0] + (my_i2-my_i1-1)*dh[0];
		xd[1] = x0[1] + (my_j2-my_j1-1)*dh[1];
		xd[2] = x0[2] + (my_k2-my_k1-1)*dh[2];

		/*number of nodes*/
		ni = my_i2-my_i1;
		nj = my_j2-my_j1;
		nk = my_k2-my_k1;
	}

	/*returns rank of processors di/dj domains from us*/
	int rankFromOffset(int di, int dj, int dk)
	{
		return (mpi_k+dk)*(mpi_size_j*mpi_size_i) + (mpi_j+dj)*mpi_size_i+(mpi_i+di);
	}

	/*returns information on MPI neighbor along a face, -1 if no neighbor*/
	int getNeighborRank(int face)
	{
		bool has_neighbor;
		int target;
		switch(face)
		{
			case 0: has_neighbor=mpi_i>0;target=rankFromOffset(-1,0,0);break;
			case 1: has_neighbor=mpi_i<mpi_size_i-1;target=rankFromOffset(1,0,0);break;
			case 2: has_neighbor=mpi_j>0;target=rankFromOffset(0,-1,0);break;
			case 3: has_neighbor=mpi_j<mpi_size_j-1;target=rankFromOffset(0,1,0);break;
			case 4: has_neighbor=mpi_k>0;target=rankFromOffset(0,0,-1);break;
			case 5: has_neighbor=mpi_k<mpi_size_k-1;target=rankFromOffset(0,0,1);break;
		}

		if (has_neighbor) return target; else return -1;
	}

	/*VARIABlES*/
	double x0[3];
	double dh[3];
	double xd[3];	/*diagonal corner*/
	int ni,nj,nk;
	
	int ***object;  /*object is used to flag internal nodes by values >0*/

	Field *phi;
	Field *rhoi;
	Field3 *ef;
	
	double dt;
	double time = 0;
	int num_colls;
	
	/*sphere data*/
	double sphere_x0[3];
	double sphere_rad;

	/*threads*/
	int num_threads = 1;

	/*MPI*/
	int mpi_rank;
	int mpi_size;
	int mpi_i,mpi_j,mpi_k;
	int mpi_size_i,mpi_size_j,mpi_size_k;
	int num_nodes_global;
};

double rnd();	/*prototype for RNG*/


#endif
