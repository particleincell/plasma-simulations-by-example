#ifndef _SOLVER_H
#define _SOLVER_H

#include "World.h"

class PotentialSolver 
{
public:
	/*constructor, sets world*/
	PotentialSolver(World &world, int max_it, double tol): 
		world(world), max_solver_it(max_it), tolerance(tol) {}
	
	/*sets reference values*/
	void setReferenceValues(double phi0, double Te0, double n0) {
		this->phi0 = phi0;
		this->Te0 = Te0;
		this->n0 = n0;
	}

	/*solves potential using Gauss-Seidel*/
	bool solve();

	/*computes electric field = -gradient(phi)*/
	void computeEF();

protected:
	World &world;
	unsigned max_solver_it;	//maximum number of solver iterations
	double tolerance;		//solver tolerance
	double phi0 = 0;		//reference plasma potential
	double n0 = 1e12;		//reference electron density
	double Te0 = 1.5;		//reference electron temperature in eV
};
#endif
