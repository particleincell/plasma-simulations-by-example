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
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

using namespace std;

/*container for storing matrix coefficients*/
struct devSeptaD
{
	double *a;
	double *b;
	double *c;
	double *d;
	double *e;
	double *f;
	double *g;
	int nu;
	int ni,nj,nk;

	/*not a destructor to avoid free after copy constructor*/
	void free () {
		cudaFree(a);
		cudaFree(b);
		cudaFree(c);
		cudaFree(d);
		cudaFree(e);
		cudaFree(f);
		cudaFree(g);
	}
};

/*builds matrix for system extended by two ghost nodes along each direction*/
void PotentialSolver::initCUDA()
{   
	threads_per_block3.x = 4;
	threads_per_block3.y = 4;
	threads_per_block3.z = 4;

	num_blocks3.x = (A.ni+threads_per_block3.x-1)/threads_per_block3.x;
	num_blocks3.y = (A.nj+threads_per_block3.y-1)/threads_per_block3.y;
	num_blocks3.z = (A.nk+threads_per_block3.z-1)/threads_per_block3.z;
	int  num_blocks = num_blocks3.x*num_blocks3.y*num_blocks3.z;

	int gni = world.ni+2;
	int gnj = world.nj+2;
	int gnk = world.nk+2;
	int gnu = gni*gnj*gnk;
	fvector a(gnu);
	fvector b(gnu);
	fvector c(gnu);
	fvector d(gnu);
	fvector e(gnu);
	fvector f(gnu);
	fvector g(gnu);

	double3 dh = world.getDh();
	float idx = 1.0/dh[0];
	float idy = 1.0/dh[1];
	float idz = 1.0/dh[2];
	float idx2 = idx*idx;	/*1/(dx*dx)*/
	float idy2 = idy*idy;
	float idz2 = idz*idz;

	/*set coefficients, loop over non-ghost nodes*/
	for (int gk=1;gk<gnk-1;gk++)
		for (int gj=1;gj<gnj-1;gj++)
			for (int gi=1;gi<gni-1;gi++)
			{
				int i = gi-1;
				int j = gj-1;
				int k = gk-1;

				int u = gk*(gni*gnj)+gj*gni+gi;

				//dirichlet node?
				if (world.object_id[i][j][k]>0)
				{
					d[u] = 1;
					continue;
				}

				//Neumann boundaries
				if (i==0) {d[u]=idx;c[u]=-idx;}
				else if (i==world.ni-1) {d[u]=idx;e[u]=-idx;}
				else if (j==0) {d[u]=idy;b[u]=-idy;}
				else if (j==world.nj-1) {d[u]=idy;f[u]=-idy;}
				else if (k==0) {d[u]=idz;a[u]=-idz;}
				else if (k==world.nk-1) {d[u]=idz;g[u]=-idz;}
				else {
					//standard internal stencil
					a[u] = idz2; g[u] = idz2;
					b[u] = idy2; f[u] = idy2;
					c[u] = idx2; e[u] = idx2;
					d[u] = -2.0*(idx2+idy2+idz2);
				}
			}

		/*allocate GPU memory for coefficient arrays*/
		cudaMalloc((void**)&devA.a, sizeof(float)*gnu);
		cudaMalloc((void**)&devA.b, sizeof(float)*gnu);
		cudaMalloc((void**)&devA.c, sizeof(float)*gnu);
		cudaMalloc((void**)&devA.d, sizeof(float)*gnu);
		cudaMalloc((void**)&devA.e, sizeof(float)*gnu);
		cudaMalloc((void**)&devA.f, sizeof(float)*gnu);
		cudaMalloc((void**)&devA.g, sizeof(float)*gnu);

		/*copy coefficients*/
		cudaMemcpy(devA.a,a,sizeof(float)*gnu,cudaMemcpyHostToDevice);
		cudaMemcpy(devA.b,b,sizeof(float)*gnu,cudaMemcpyHostToDevice);
		cudaMemcpy(devA.c,c,sizeof(float)*gnu,cudaMemcpyHostToDevice);
		cudaMemcpy(devA.d,d,sizeof(float)*gnu,cudaMemcpyHostToDevice);
		cudaMemcpy(devA.e,e,sizeof(float)*gnu,cudaMemcpyHostToDevice);
		cudaMemcpy(devA.f,f,sizeof(float)*gnu,cudaMemcpyHostToDevice);
		cudaMemcpy(devA.g,g,sizeof(float)*gnu,cudaMemcpyHostToDevice);

		/*copy matrix to

		/*allocate memory for potential and charge density*/
		cudaMalloc((void**)&dev_phi, sizeof(float)*gnu);
		cudaMalloc((void**)&dev_b, sizeof(float)*gnu);
		cudaMalloc((void**)&dev_res, num_blocks*sizeof(float));

		/*allocate CPU memory for res*/
		cudaHostAlloc((void**)&host_res,num_blocks*sizeof(float),cudaHostAllocDefault);
}


/*Gauss-Seidel Poisson solver*/

/*electron reference parameters*/
__constant__ double dev_n0;
__constant__ double dev_phi0;
__constant__ double dev_kTe0;

/*cuda kernel*/
__global__ void cudaGSupdate (devSeptaD *A, double *phi, double *b, char *object)
{
	int i = blockIdx.x*blockDim.x+threadIdx.x;
	int j = blockIdx.y*blockDim.y+threadIdx.y;
	int k = blockIdx.z*blockDim.z+threadIdx.z;
	int ni = A->ni;
	int nj = A->nj;
	int nk = A->nk;

	/*compute index*/
	int u = k*ni*nj+j*ni+i;
	
	if (i>0 && i<ni-1 &&
		j>0 && j<nj-1 &&
		k>0 && k<nk-1)
	{
		double rhoe = 0;
		
		/*open node*/	
		if (object[u]==0)
			rhoe = (QE*dev_n0*exp((phi[u] - dev_phi0)/dev_kTe0))/EPS_0;

		double g = ((b[u] + rhoe) -
					A->a[u]*phi[u-ni*nj] -
					A->b[u]*phi[u-ni] -
					A->c[u]*phi[u-1] -
					A->e[u]*phi[u+1] -
					A->f[u]*phi[u+ni] -
					A->g[u]*phi[u+ni*nj])/A->d[u];

		/*SOR not converging with Jacobi*/
		phi [u] = g;		
	}
}

/*convergence check*/
__global__ void cudaGSresidue(double *res, devSeptaD *A, double *phi, double *b, char *object)
{
	__shared__ float my_res[1024];		/*1024 is max threads per block*/
	double R = 0;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int tz = threadIdx.z;
	int i = blockIdx.x*blockDim.x+tx;
	int j = blockIdx.y*blockDim.y+ty;
	int k = blockIdx.z*blockDim.z+tz;
	int ni = A->ni;
	int nj = A->nj;
	int nk = A->nk;

	/*compute index*/
	int u = k*ni*nj+j*ni+i;
	
	if (i>0 && i<ni-1 &&
		j>0 && j<nj-1 &&
		k>0 && k<nk-1)
	{
		double rhoe = 0;
        
        if (object[u]==0)
            rhoe = (QE*dev_n0*exp((phi[u] - dev_phi0)/dev_kTe0))/EPS_0;
						
	   R = (b[u] + rhoe) -
            A->a[u]*phi[u-ni*nj] - 
            A->b[u]*phi[u-ni] - 
            A->c[u]*phi[u-1] - 
            A->d[u]*phi[u] -
            A->e[u]*phi[u+1] - 
            A->f[u]*phi[u+ni] - 
            A->g[u]*phi[u+ni*nj];
	}

	my_res[tz*blockDim.x*blockDim.y+ty*blockDim.x+tx] = R*R;

	/*wait for all threads from block to finish*/
	__syncthreads();

	/*if this is "root", sum up, slow way*/
	if (tx==0 && ty==0 && tz==0)
	{
		double sum = 0;
		for (int i=0;i<blockDim.x*blockDim.y*blockDim.z;i++)
		{
			sum+=my_res[i];
		}
		/*save in global memory*/
		res[blockIdx.z*gridDim.x*gridDim.y+
			blockIdx.y*gridDim.x+
			blockIdx.x] = sum;
	}            
}

/*updated version that leaves potential on the GPU*/
bool PotentialSolver::solveGSCUDA()
{
	bool converged = false;
	double L2;

	/*copy potential on the first time*/
	if (first_time)
	{		
		deflate(phi,world.phi.data);
		CUDA_ERROR(cudaMemcpy(dev_phi,phi,A.nu*sizeof(double),cudaMemcpyHostToDevice));
		first_time = false;
	}

	/*compute number of blocks for residue checking*/
	int num_blocks = num_blocks3.x*num_blocks3.y*num_blocks3.z;
	 
    /*set RHS to zero on boundary nodes (zero electric field) 
    and to existing potential on fixed nodes        */
	deflate(b,world.rhoi->data);
	
	for (int u=0;u<A.nu;u++)
    {
		if (object[u]<0) b[u] = 0;				/*neumann boundary*/
        else if (object[u]>0) b[u] = phi[u];	/*dirichlet boundary*/
        else b[u] *= -1.0/EPS_0;                /*open node*/
    }
        	
	/*now copy data*/
	CUDA_ERROR(cudaMemcpy(dev_b,b,A.nu*sizeof(double),cudaMemcpyHostToDevice));
		
	/*solve potential*/
	
	int solver_it;

	//int max_it=25;
	for (solver_it=0;solver_it<max_it;solver_it++)
	{
		/*launch threads*/
		cudaGSupdate<<<num_blocks3,threads_per_block3>>>(dev_devA,dev_phi,dev_b);
	
		if (solver_it%25==0)
		{
			cudaGSresidue<<<num_blocks3,threads_per_block3>>>(dev_res,dev_devA,dev_phi,dev_b,dev_object);
			cudaMemcpy(res_pinned,dev_res,num_blocks*sizeof(double),cudaMemcpyDeviceToHost);

			double sum=0;
			for (int i=0;i<num_blocks;i++) sum+=res_pinned[i];
		    L2 = sqrt(sum/(A.nu));
			if (L2<tol) {converged=true;break;}
		}
	}

	/*we leave potential on the GPU so don't need to copy back*/
	if (!converged) cerr<<"cudaGS failed to converge, L2 = "<<L2<<endl;
	
    return converged;
}

/********* ADDITIONAL CPU CODE TO SUPPORT CUDA ********************/
/*memory cleanup*/
PotentialSolver::~PotentialSolver()
{
	if (solver_type==GSCUDA)
	{
		devA.free();
		delete(devA);
		CUDA_ERROR(cudaFree(dev_devA));
	
		CUDA_ERROR(cudaFree(dev_phi));
		CUDA_ERROR(cudaFree(dev_b));
		CUDA_ERROR(cudaFree(dev_object));
		CUDA_ERROR(cudaFree(dev_res));

		/*also free host memory*/
		CUDA_ERROR(cudaFreeHost(res_pinned));
		delete[] phi;
		delete[] b;	
		CUDA_ERROR(cudaFree(world.dev_ef3));
	}

}

/*updates phi on the CPU*/
void PotentialSolver::updateHostPhi()
{
	if (solver_type!=GSCUDA) return;

	cudaMemcpy(phi,dev_phi,A.nu*sizeof(double),cudaMemcpyDeviceToHost);
	inflate(phi, world.phi.data);
}
