/*solves poisson equation with Boltzmann electrons using the Gauss-Seidel scheme*/
#ifndef _SOLVER_H
#define _SOLVER_H

#include "World.h"
#include <math.h>
#include <string.h>
#include <iostream>
#include <thread>
#include <condition_variable>

/*septa-diagonal matrix*/
class SeptaD
{
public: 
    SeptaD(World &world);
	~SeptaD() {
		delete[] a;
		delete[] b;
		delete[] c;
		delete[] d;
		delete[] e;
		delete[] f;
		delete[] g;
	}
	/*vectors storing matrix data*/
	double  *a,*b,*c,*d,*e,*f,*g;
	int nu;
  	int ni,nj,nk;
};

/*from http://stackoverflow.com/questions/24465533/implementing-boostbarrier-in-c11*/
class Barrier {
public:
    explicit Barrier(std::size_t iCount) :
      mThreshold(iCount),
      mCount(iCount),
      mGeneration(0) {
    }

    void Wait() {
        auto lGen = mGeneration;
        std::unique_lock<std::mutex> lLock{mMutex};
        if (!--mCount) {
            mGeneration++;
            mCount = mThreshold;
            mCond.notify_all();
        } else {
            mCond.wait(lLock, [this, lGen] { return lGen != mGeneration; });
        }
    }

private:
    std::mutex mMutex;
    std::condition_variable mCond;
    std::size_t mThreshold;
    std::size_t mCount;
    std::size_t mGeneration;
};


enum SolverType {GS, GSMT, QN};

/*potential solver class*/
class PotentialSolver 
{
public:
	/*constructor, sets world*/
	PotentialSolver(World &world, SolverType type): 
		world(world), A(world), solver_type(type),
		barrier(world.num_threads){
		ni = A.ni;
		nj = A.nj;
		nk = A.nk;
		object = new int[A.nu];
		max_it = 2500;
		tol = 1e-1;
		my_sum = new double[world.num_threads];

		buildMatrix();
	}

	~PotentialSolver() {delete[] object;delete[] my_sum;}
	
	/*builds the "A" matrix for linear potential solver*/
    void buildMatrix();

	/*sets reference values for Boltzmann electron relationship*/
	void setElectronRef(double phi0, double n0, double kTe0)
	{
		this->phi0 = phi0;
		this->n0 = n0;
		this->kTe0 = kTe0;
	}
	
	bool solvePotential()
    {
        switch(solver_type)
        {
            case GS: return solveGS();
            case GSMT: return solveGSMT();
            case QN: return solveQN();
            default: return false;
        }
    }
    
 	/*computes electric field = -gradient(phi)*/
	void computeEF();

	/** convert 3D data to a 1D vector*/
    void deflate(double *r, double ***data)
    {
    	/*set all to zero*/
    	memset(r,0,sizeof(double)*A.nu);

        for (int i=1;i<ni-1;i++)
                for (int j=1;j<nj-1;j++)
                    for (int k=1;k<nk-1;k++)
                         r[ijkn(i,j,k)] = data[i-1][j-1][k-1];
	}
	
    /** convert 1D vector to 2D array*/
    void inflate(double ***data3d, double *data1d)
    {
            for (int i=1;i<ni-1;i++)
                    for (int j=1;j<nj-1;j++)
                        for (int k=1;k<nk-1;k++)
                            data3d[i-1][j-1][k-1] = data1d[ijkn(i,j,k)];
    }

	int ijkn(int i, int j, int k) {return k*ni*nj+j*ni+i;}

	static void gsMTKernel(int thread_id, PotentialSolver *solver, double *phi, double *b);

protected:
	/*computes potential using quasineutral boltzmann model*/
    bool solveQN();

	/*solves potential using Gauss-Seidel*/
	bool solveGS();
	bool solveGSMT();

	/*ghost nodes update for MPI runs*/
	void updateGhosts(double *phi);

	int max_it;	/*number of solver iterations*/
    double tol;	/*solver tolerance*/

	World &world;
	double phi0,n0,kTe0;		/*reference values for electron model*/
	SeptaD A;
	SolverType solver_type;
	int *object;
	int ni,nj,nk;

	/*variables used by threads*/
	Barrier barrier;
	double *my_sum;

	/*mpi*/
	int rank[6];	/*ranks of MPI neighbors*/
};




#endif
