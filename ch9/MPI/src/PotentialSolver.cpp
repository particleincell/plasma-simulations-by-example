/*solves poisson equation with Boltzmann electrons using the Gauss-Seidel scheme*/

#include "PotentialSolver.h"
#include "Field.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <thread>
#include <stdexcept>
#include <mpi.h>
#include "World.h"

using namespace std;

/*vector math helper functions*/
namespace vec
{
	/*sets vector to zero*/
	void clear(double v[], int nu)
	{
		for (int u=0;u<nu;u++) v[u]=0;
	}

	/*allocates and clears an array*/
	double* alloc(int nu)
	{
		double *r =new double[nu];
		clear(r,nu);
		return r;
	}
};


/*SeptaD constructor*/
SeptaD::SeptaD(World &world)
{
	/*add ghost cells*/
	ni = world.ni+2;
	nj = world.nj+2;
	nk = world.nk+2;

	nu = ni*nj*nk;

    a = vec::alloc(nu);
    b = vec::alloc(nu);
    c = vec::alloc(nu);
    d = vec::alloc(nu);
    e = vec::alloc(nu);
    f = vec::alloc(nu);
    g = vec::alloc(nu);
}

/*builds matrix for system extended by two ghost nodes along each direction*/
void PotentialSolver::buildMatrix()
{   
	double idx = 1.0/world.dh[0];
	double idy = 1.0/world.dh[1];
	double idz = 1.0/world.dh[2];
    double idx2 = 1.0/(world.dh[0]*world.dh[0]);	/*1/(dx*dx)*/
	double idy2 = 1.0/(world.dh[1]*world.dh[1]);
	double idz2 = 1.0/(world.dh[2]*world.dh[2]);

	/*set ranks of MPI neighbors*/
	for (int i=0;i<6;i++) rank[i] = world.getNeighborRank(i);

	/*initialize all types to "fixed" to mark ghost nodes,
	 * will be overwritten on internal [1:ni-1]..*/
	for (int u=0;u<A.nu;u++) {A.d[u]=1; object[u]=1;}

	/*process internal, non ghost nodes*/
    for (int i=1;i<ni-1;i++)
        for (int j=1;j<nj-1;j++)
            for (int k=1;k<nk-1;k++)
            {
            	/*index in our shifted system*/
                int u = ijkn(i,j,k);
                
                /*object using global indexing*/
				object[u] = world.object[i-1][j-1][k-1];
                if (object[u]>0)
                {
                    A.d[u]=1;
                    continue;
                }

                /*set neuman boundaries on non-processor boudaries*/
                /*note shift in index from i==0 to i==1 and so on*/
                if (i==1 && rank[0]<0) {A.d[u]=idx;A.e[u]=-idx;object[u]=-1;}
                else if (i==ni-2 && rank[1]<0) {A.d[u]=idx;A.c[u]=-idx;object[u]=-1;}
                else if (j==1 && rank[2]<0) {A.d[u]=idy;A.f[u]=-idy;object[u]=-1;}
                else if (j==nj-2 && rank[3]<0) {A.d[u]=idy;A.b[u]=-idy;object[u]=-1;}
                else if (k==1 && rank[4]<0) {A.d[u]=idz;A.g[u]=-idz;object[u]=-1;}
				else if (k==nk-2 && rank[5]<0) {A.d[u]=idz;A.a[u]=-idz;object[u]=-1;}
                else	/*standard central difference*/
				{
                    A.a[u] = idz2;
                    A.b[u] = idy2;
                    A.c[u] = idx2;
                    A.d[u] = -2.0*(idx2+idy2+idz2);
                    A.e[u] = idx2;
                    A.f[u] = idy2;
                    A.g[u] = idz2;
                }
				
					
            }
}

/*solves potential using Gauss-Seidel*/
bool PotentialSolver::solveGS()
{
 	double *phi = new double[A.nu];
	double *b = new double[A.nu];

	auto clock_start = chrono::high_resolution_clock::now();

	deflate(phi,world.phi->data);
	deflate(b,world.rhoi->data);
     
    /*set RHS to zero on boundary nodes (zero electric field) 
    and to existing potential on fixed nodes        */
    for (int u=0;u<ni*nj*nk;u++)
    {
		if (object[u]<0) b[u] = 0;				/*neumann boundary*/
        else if (object[u]>0) b[u] = phi[u];	/*dirichlet boundary*/
        else b[u] = -b[u]/EPS_0;                /*open node*/
    }
        
	double L2=0;
	bool converged = false;
		
	/*solve potential*/
	for (int solver_it=0;solver_it<max_it;solver_it++)
	{
		/*update ghost nodes*/
		updateGhosts(phi);

		for (int u=0;u<ni*nj*nk;u++)
        {
			double rhoe = 0;

			/*added for MPI, otherwise ghosts will get replaced by b[u]*/
			if (object[u]>0) continue;

			if (object[u]==0)
				rhoe = (QE*n0*exp((phi[u] - phi0)/kTe0))/EPS_0;

			double phi_new = ((b[u] + rhoe) -
					(A.a[u]!=0?A.a[u]*phi[u-ni*nj]:0) -
					(A.b[u]!=0?A.b[u]*phi[u-ni]:0) -
					(A.c[u]!=0?A.c[u]*phi[u-1]:0) -
					(A.e[u]!=0?A.e[u]*phi[u+1]:0) -
					(A.f[u]!=0?A.f[u]*phi[u+ni]:0) -
					(A.g[u]!=0?A.g[u]*phi[u+ni*nj]:0))/A.d[u];

            /*SOR*/
            phi[u] = phi[u] + 1.4*(phi_new-phi[u]);
        }

		/*check for convergence*/
		if (solver_it%25==0)
		{
			double my_sum = 0;

            for (int u=0;u<ni*nj*nk;u++)
            {
                double rhoe = 0;
                
            	if (object[u]>0) continue;

                if (object[u]==0)
                    rhoe = (QE*n0*exp((phi[u] - phi0)/kTe0))/EPS_0;
						
                double R = (b[u] + rhoe) -
                (A.a[u]!=0?A.a[u]*phi[u-ni*nj]:0) - 
                (A.b[u]!=0?A.b[u]*phi[u-ni]:0) - 
                (A.c[u]!=0?A.c[u]*phi[u-1]:0) - 
                (A.d[u]!=0?A.d[u]*phi[u]:0) -
                (A.e[u]!=0?A.e[u]*phi[u+1]:0) - 
                (A.f[u]!=0?A.f[u]*phi[u+ni]:0) - 
                (A.g[u]!=0?A.g[u]*phi[u+ni*nj]:0);
                        
				my_sum+=R*R;
			}
			
            /*compute global sum*/
            double glob_sum;
            MPI_Allreduce(&my_sum,&glob_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

			L2 = sqrt(glob_sum/world.num_nodes_global);
			if (L2<tol) {converged=true;break;}
		}
	}
		
    inflate(world.phi->data,phi);

	delete[] phi;
	delete[] b;
        
	if (!converged && world.mpi_rank==0) cerr<<"GS failed to converge, L2 = "<<L2<<endl;

	/*show timing data?*/
	if (false)
	{
		auto clock_end = chrono::high_resolution_clock::now();
		std::chrono::duration<float> delta = clock_end-clock_start;

		if (world.mpi_rank==0)
			cout << "GS took "<<delta.count()<< "s\n";
	}
    return converged;
}

/*performs GS update on on a chunk*/
void PotentialSolver::gsMTKernel(int thread_id, PotentialSolver *solver, 
double *phi_in, double *b)
{
	/*grab data*/
	double n0 = solver->n0;
	double phi0 = solver->phi0;
	double kTe0 = solver->kTe0;
	int *object = solver->object;
	SeptaD &A = solver->A;
	int ni = solver->ni;
	int nj = solver->nj;
	int nk = solver->nk;
	double *phi = phi_in;
	
	/*figure out our own extend using thread id*/
	int chunk = A.nu/solver->world.num_threads + 1;
	int u1 = thread_id*chunk;
	int u2 = (thread_id+1)*chunk;
	if (u2>A.nu) u2=A.nu;
	solver->my_sum[thread_id]=1e66;		/*init to some large value*/

	for (int solver_it = 0;solver_it<solver->max_it;solver_it++)
	{
		/*compute new potential*/
		for (int u=u1;u<u2;u++)
		{
			double rhoe = 0;

			if (object[u]==0)
				rhoe = (QE*n0*exp((phi[u] - phi0)/kTe0))/EPS_0;

			double phi_new = ((b[u] + rhoe) -
					(A.a[u]!=0?A.a[u]*phi[u-ni*nj]:0) -
					(A.b[u]!=0?A.b[u]*phi[u-ni]:0) -
					(A.c[u]!=0?A.c[u]*phi[u-1]:0) -
					(A.e[u]!=0?A.e[u]*phi[u+1]:0) -
					(A.f[u]!=0?A.f[u]*phi[u+ni]:0) -
					(A.g[u]!=0?A.g[u]*phi[u+ni*nj]:0))/A.d[u];

	        /*SOR, check slides on this*/
            phi[u] = phi[u] + 1.4*(phi_new-phi[u]);
   		}

		/*thread barrier, wait for all threads to finish*/
		solver->barrier.Wait();

		/*do this the following iteration after computing residue*/
		if (solver_it%25==1)
		{
			/*convergence check, figure out global L2*/
			double global_sum = 0;
			for (int i=0;i<solver->world.num_threads;i++)
				global_sum += solver->my_sum[i];

			double L2 = sqrt(global_sum/(ni*nj*nk));
			//if (thread_id==0) cout<<L2<<endl;
			if (L2<solver->tol)	break;
		}

		/*check for convergence*/
		if (solver_it%25==0)
		{
			double sum = 0;

			/*note that each thread is computing residue only for its chunk*/
			for (int u=u1;u<u2;u++)
			{
				double rhoe = 0;

				if (object[u]==0)
					rhoe = (QE*n0*exp((phi[u] - phi0)/kTe0))/EPS_0;

				double R = (b[u] + rhoe) -
				(A.a[u]!=0?A.a[u]*phi[u-ni*nj]:0) -
				(A.b[u]!=0?A.b[u]*phi[u-ni]:0) -
				(A.c[u]!=0?A.c[u]*phi[u-1]:0) -
				(A.d[u]!=0?A.d[u]*phi[u]:0) -
				(A.e[u]!=0?A.e[u]*phi[u+1]:0) -
				(A.f[u]!=0?A.f[u]*phi[u+ni]:0) -
				(A.g[u]!=0?A.g[u]*phi[u+ni*nj]:0);

				sum+=R*R;
			}

			solver->my_sum[thread_id] = sum;

		}
	} /*for loop*/
}

/*solves potential using multi-threaded GS
 * note this version is not compatible with MPI (yet!)
 * */
bool PotentialSolver::solveGSMT()
{
 	double *phi = new double[A.nu];
	double *b = new double[A.nu];
	bool converged = false;

	if (world.mpi_size>1) {throw runtime_error("solveGSMT not yet compatible with MPI!");}

	deflate(phi,world.phi->data);
	deflate(b,world.rhoi->data);

    /*set RHS to zero on boundary nodes (zero electric field)
    and to existing potential on fixed nodes        */
    for (int u=0;u<ni*nj*nk;u++)
    {
		if (object[u]<0) b[u] = 0;				/*neumann boundary*/
        else if (object[u]>0) b[u] = phi[u];	/*dirichlet boundary*/
        else b[u] = -b[u]/EPS_0;                /*open node*/
    }

	thread *updaters = new thread[world.num_threads];
	for (int i=0;i<world.num_threads;i++)
	{
		/*launch new thread*/
		updaters[i] = thread(gsMTKernel,i,this, phi, b);
	}

	/*wait for threads to finish*/
	for (int i=0;i<world.num_threads;i++)
		updaters[i].join();

	/*figure out global L2*/
	double sum=0;
	for (int i=0;i<world.num_threads;i++)
		sum += my_sum[i];

	double L2 = sqrt(sum/(ni*nj*nk));
	if (L2<tol)	converged = true;
	if (!converged) cerr<<"gsMT failed to converge, L2 = "<<L2<<endl;

    inflate(world.phi->data,phi);
    delete[] phi;
	delete[] b;
	delete[] updaters;

	return converged;
}

/*quasi-neutral potential solver*/
bool PotentialSolver::solveQN()
{
	double ***phi = world.phi->data;
    double ***rhoi = world.rhoi->data;
    double rho0 = n0*QE;
    double phi_min = phi0 + kTe0*log(1e-6);
        
    for (int i=0;i<ni;i++)
        for (int j=0;j<nj;j++)
            for (int k=0;k<nk;k++)
            {
                if (world.object[i][j][k]>0) continue; /*skip solids*/
                    
                if (rhoi[i][j][k]>0)
                    phi[i][j][k] = phi0 + kTe0*log(rhoi[i][j][k]/rho0);
                else 
                    phi[i][j][k] = phi_min;
            }
    return true;
}

/*computes electric field = -gradient(phi)*/
void PotentialSolver::computeEF()
{
	int ni = world.ni;
	int nj = world.nj;
	int nk = world.nk;
	double dx = world.dh[0];
	double dy = world.dh[1];
	double dz = world.dh[2];
		
	double ***phi = world.phi->data;
	double ***ex = world.ef->f[0]->data;
	double ***ey = world.ef->f[1]->data;
	double ***ez = world.ef->f[2]->data;
	
	for (int i=0;i<ni;i++)
		for (int j=0;j<nj;j++)
			for (int k=0;k<nk;k++)
			{
				/*x component*/
				if (i==0) 
					ex[i][j][k] = -(phi[i+1][j][k]-phi[i][j][k])/dx;	/*forward*/
				else if (i==ni-1) 
					ex[i][j][k] = -(phi[i][j][k]-phi[i-1][j][k])/dx;	/*backward*/
				else
					ex[i][j][k] = -(phi[i+1][j][k] - phi[i-1][j][k])/(2*dx);	/*central*/

				/*y component*/
				if (j==0) 
					ey[i][j][k] = -(phi[i][j+1][k]-phi[i][j][k])/dy;
				else if (j==nj-1) 
					ey[i][j][k] = -(phi[i][j][k]-phi[i][j-1][k])/dy;
				else
					ey[i][j][k] = -(phi[i][j+1][k] - phi[i][j-1][k])/(2*dy);

				/*z component*/
				if (k==0) 
					ez[i][j][k] = -(phi[i][j][k+1]-phi[i][j][k])/dz;
				else if (k==nk-1) 
					ez[i][j][k] = -(phi[i][j][k]-phi[i][j][k-1])/dz;
				else
					ez[i][j][k] = -(phi[i][j][k+1] - phi[i][j][k-1])/(2*dz);
			}
}

/*updates ghost nodes, similar to Field::updateBoundaries*/
void PotentialSolver::updateGhosts(double *phi)
{
	const int TAG_DATA = 20;

	for (int face=0;face<6;face++)
	{
		int target = rank[face];

		/*transfer data*/
		if (target>=0)
		{
			int size;
			if (face/2==0) size=(nj-2)*(nk-2);
			else if (face/2==1) size=(ni-2)*(nk-2);
			else if (face/2==2) size=(ni-2)*(nj-2);

			/*copy values along left boundary to send buffer*/
			double *send_buffer = new double[size];

			/*pack data, there is probably a more clever way of doing this*/
			double *s=send_buffer;	/*pointer to start*/
			switch (face)
			{
				case 0:  for (int j=1;j<nj-1;j++) for (int k=1;k<nk-1;k++) (*s++)=phi[ijkn(2,j,k)];break;
				case 1:  for (int j=1;j<nj-1;j++) for (int k=1;k<nk-1;k++) (*s++)=phi[ijkn(ni-3,j,k)];break;
				case 2:  for (int i=1;i<ni-1;i++) for (int k=1;k<nk-1;k++) (*s++)=phi[ijkn(i,2,k)];break;
				case 3:  for (int i=1;i<ni-1;i++) for (int k=1;k<nk-1;k++) (*s++)=phi[ijkn(i,nj-3,k)];break;
				case 4:  for (int i=1;i<ni-1;i++) for (int j=1;j<nj-1;j++) (*s++)=phi[ijkn(i,j,2)];break;
				case 5:  for (int i=1;i<ni-1;i++) for (int j=1;j<nj-1;j++) (*s++)=phi[ijkn(i,j,nk-3)];break;
			}

			double *recv_buffer = new double[size];

			/*transfer values*/
			MPI_Sendrecv(send_buffer, size, MPI_DOUBLE, target, TAG_DATA,
						 recv_buffer, size, MPI_DOUBLE, target, TAG_DATA,
						 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			/*update our ghost cells, note difference in indexes from send pack*/
			double *r=recv_buffer;	/*pointer to start*/
			switch (face)
			{
				case 0:  for (int j=1;j<nj-1;j++) for (int k=1;k<nk-1;k++) phi[ijkn(0,j,k)]=(*r++);break;
				case 1:  for (int j=1;j<nj-1;j++) for (int k=1;k<nk-1;k++) phi[ijkn(ni-1,j,k)]=(*r++);break;
				case 2:  for (int i=1;i<ni-1;i++) for (int k=1;k<nk-1;k++) phi[ijkn(i,0,k)]=(*r++);break;
				case 3:  for (int i=1;i<ni-1;i++) for (int k=1;k<nk-1;k++) phi[ijkn(i,nj-1,k)]=(*r++);break;
				case 4:  for (int i=1;i<ni-1;i++) for (int j=1;j<nj-1;j++) phi[ijkn(i,j,0)]=(*r++);break;
				case 5:  for (int i=1;i<ni-1;i++) for (int j=1;j<nj-1;j++) phi[ijkn(i,j,nk-1)]=(*r++);break;
			}

			delete[] send_buffer;
			delete[] recv_buffer;
		}
	} /*for face*/
};
