#ifndef _SOLVER_H
#define _SOLVER_H

#include <assert.h>
#include "World.h"

enum SolverType {GS, QN};

class PotentialSolver
{
public:
	/*constructor*/
	PotentialSolver(World &world, SolverType type, int max_it, double tol):
		world(world), solver_type(type),
		max_solver_it(max_it), tolerance(tol) {	}

	/*sets reference values*/
	void setReferenceValues(double phi0, double Te0, double n0) {
		this->phi0 = phi0;
		this->Te0 = Te0;
		this->n0 = n0;
	}

	/*computes electric field = -gradient(phi)*/
	void computeEF();

	/*builds the "A" matrix for linear potential solver*/
    void buildMatrix();

    //calls the appropriate potential solver
    bool solve()
    {
    	switch(solver_type)
    	{
    		case GS: return solveGS();
    		case QN: return solveQN();
    		default: return false;
    	}
    }
protected:
	World &world;
    SolverType solver_type;

    enum NodeType {REG,NEUMANN,DIRICHLET};
    std::vector<NodeType> node_type;	//flag for different node types

	unsigned max_solver_it;	//maximum number of solver iterations
	double tolerance;		//solver tolerance
	double phi0 = 0;		//reference plasma potential
	double n0 = 1e12;		//reference electron density
	double Te0 = 1.5;		//reference electron temperature in eV

	/*computes potential using quasineutral boltzmann model*/
	bool solveQN();

	/*solves non-linear potential using Gauss-Seidel*/
	bool solveGS();
};
#endif
