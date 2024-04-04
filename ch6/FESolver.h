#ifndef _FESOLVER_H
#define _FESOLVER_H

#include "World.h"

//useful utility functions
namespace utils {
	/*computes determinant of a 3x3 matrix*/
	double det3(double (*M)[3]);

	/*computes determinant of a 4x4 matrix*/
	double det4(double (*M)[4]);
};


//structure to hold data for a single row
template <int S>
struct Row {
	Row() {for (int i=0;i<S;i++) {a[i]=0;col[i]=-1;}}
	void operator= (const Row &o) {for (int i=0;i<S;i++) {a[i] = o.a[i];col[i]=o.col[i];}}
	double a[S];		//coefficients
	int col[S];
};


/*solver class*/
class FESolver
{
public:
	
	FESolver(World &world, int max_it, double tol); /*constructor, initialized data structures*/
	~FESolver();	/*destructor, frees memory*/

	/*sets reference values*/
	void setReferenceValues(double phi0, double Te0, double n0) {
		this->phi0 = phi0;
		this->Te0 = Te0;
		this->n0 = n0;
	}

	void solve();
	void computeEF();

	/*evaluates ef in a cell. Since constant field in cell, just copy*/
	double3 evalEf(LCord &lc) {return ef[lc.cell_id];}

protected:
	void computeNX();
	void startAssembly();	/*clears K and F*/
	void preAssembly();
	void addKe(int e, double ke[4][4]);	/*adds contributions from element stiffness matrix*/
	void addFe(dvector &F, int e, double fe[4]); /*adds contributions from element force vector*/

	double evalNa(int a, double xi, double eta, double zeta);
	double3 getNax(int e, int a, double xi, double eta, double zeta);

	void inverse(double M[3][3], double V[3][3]);
	void buildF1Vector(dvector &ion_den);
	void solveNonLinear(dvector &d);
	void solveLinear(double **K, dvector &d, dvector &F);	/*solves Kd=F for d*/

	void buildJmatrix();
	
	/*helper functions for matrix math, y=A*x */
	void matVecMultiply(dvector &y, double **A, dvector &x, int nu)
	{
		for (int i=0;i<nu;i++)
		{
			y[i] = 0;
			for (int j=0;j<nu;j++)
				y[i] += A[i][j]*x[j];
		}
	}

	/*computes y=v1-v2*/
	void vecVecSubtract(dvector &y, dvector &v1, dvector &v2, int nu)
	{
		for (int i=0;i<nu;i++)
				y[i] = v1[i]-v2[i];
	}

	World &world;
	VolumeMesh &vm;
	int n_nodes;
	int n_elements;	/*save this so we can properly deallocate LM*/

	double **K;		/*global stiffness matrix, should use a sparse matrix*/
	double **J;		/*Jacobian matrix*/
	int **LM;		/*LM[e][a] location matrix */
	double ***NX;	/*NX[e][a] is a dNa/dx [3] vector*/

	dvector F0;		/*"fh" and "fg" parts of the global force vector*/
	dvector F1;		/*"ff" part of the force vector*/
	ivector ID;		/*ID[n]=A*/
	int neq;		/*number of unknowns/equations*/

	/*solution*/
	dvector d;		/*d[neq] is the solution on the unknown nodes*/
	dvector g;		/*g[n] essential boundaries*/
	dvector uh;		/*uh[n] solution on nodes, union of d and g*/

	dvector3 ef;	/*ef[e] is the electric field in cell e*/

	dvector detJ; /*determinant of the jacobian x_xi*/

	/*quadrature points*/
	double l[2];
	double W[2];
	int n_int;

	unsigned max_solver_it;	//maximum number of solver iterations
	double tolerance;		//solver tolerance
	double phi0 = 0;		//reference plasma potential
	double n0 = 1e12;		//reference electron density
	double Te0 = 1.5;		//reference electron temperature in eV
};


#endif
