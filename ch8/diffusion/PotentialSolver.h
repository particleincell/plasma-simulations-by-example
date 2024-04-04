#ifndef _SOLVER_H
#define _SOLVER_H

#include <assert.h>
#include "World.h"

//structure to hold data for a single row
template <int S>
struct Row {
	Row() {for (int i=0;i<S;i++) {a[i]=0;col[i]=-1;}}
	void operator= (const Row &o) {for (int i=0;i<S;i++) {a[i] = o.a[i];col[i]=o.col[i];}}
	double a[S];		//coefficients
	int col[S];
};

/*matrix with up to seven non zero diagonals*/
class Matrix
{
public:
    Matrix(int nr):nu{nr} {rows=new Row<nvals>[nr];}
    Matrix(const Matrix &o):Matrix(o.nu) {
    	for (int r=0;r<nu;r++) rows[r] = o.rows[r];
    };	//copy constructor
    ~Matrix() {if (rows) delete[] rows;}
	dvector operator*(dvector &v);	//matrix-vector multiplication

	double& operator() (int r, int c); //reference to A[r,c] value in a full matrix
	void clearRow(int r) {rows[r]=Row<nvals>();} //reinitializes a row
	Matrix diagSubtract(dvector &P);	//subtracts a vector from the diagonal
	Matrix invDiagonal();		//returns a matrix containing inverse of our diagonal
	double multRow(int r, dvector &x);	//multiplies row r with vector x

	static constexpr int nvals = 7;	//maximum 7 non-zero values
	const int nu;			//number of rows (unknowns)

protected:
	Row<nvals> *rows;	//row data
};

enum SolverType {GS, PCG, QN};

class PotentialSolver
{
public:
	/*constructor*/
	PotentialSolver(World &world, SolverType type, int max_it, double tol):
		world(world), solver_type(type), A(world.ni*world.nj*world.nk),
		max_solver_it(max_it), tolerance(tol) {
			buildMatrix();
		}

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
    		case PCG: return solveNRPCG();
    		case QN: return solveQN();
    		default: return false;
    	}
    }
protected:
	World &world;
    SolverType solver_type;
    Matrix A;				//system matrix for the linear equation

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

	/*linear PCG solver for Ax=b system*/
	bool solvePCGLinear(Matrix &A, dvector &x, dvector &b);

	/*linear GS solver for Ax=b system*/
	bool solveGSLinear(Matrix &A, dvector &x, dvector &b);

	/*Newton Raphson solver for a nonlinear system, uses PCG for the linear solve*/
	bool solveNRPCG();
};
#endif
