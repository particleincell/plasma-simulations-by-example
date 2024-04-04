#ifndef _MAGNETIC_SOLVER_H
#define _MAGNETIC_SOLVER_H

#include <assert.h>
#include "World.h"
#include "PotentialSolver.h"

class MagneticSolver
{
public:
	/*constructor*/
	MagneticSolver(World &world, int max_it, double tol):
		world(world), A(world.ni*world.nj*world.nk),
		max_solver_it(max_it), tolerance(tol) {
			buildMatrix();
		}

	/*builds the "A" matrix for linear potential solver*/
    void buildMatrix();

    //solves magnetic potential
    bool solve();

    //sets the analytical solution per Jackson
    void analyticalSol();

	//computes a flattened version of a divergence field
	dvector computeDivergence(Field3 &M);

protected:
	World &world;
    Matrix A;				//system matrix for the linear equation
    dvector phi_m;
    int sphere_orig_u;

    enum NodeType {REG,NEUMANN,DIRICHLET};
    std::vector<NodeType> node_type;	//flag for different node types

	unsigned max_solver_it;	//maximum number of solver iterations
	double tolerance;		//solver tolerance
};
#endif
