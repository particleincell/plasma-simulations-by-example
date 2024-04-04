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

	static constexpr int nvals = 5;	//maximum 7 non-zero values
	const int nu;			//number of rows (unknowns)

protected:
	Row<nvals> *rows;	//row data
};


/*vector math helper functions*/
namespace vec
{
	/*returns sum of v1[i]*v2[i]*/
	double dot(dvector v1, dvector v2);
	/*returns l2 norm*/
	double norm(dvector v);

	/** converts 3D field to a 1D vector*/
	dvector deflate(Field &f3);

	/** converts 1D vector to 3D field*/
	void inflate(dvector &d1, Field& f3);

	/** converts 1D vector to 3D field*/
	Field inflate(dvector &d1,int ni, int nj);
};


class EMSolver
{
public:
	/*constructor*/
	EMSolver(World &world, int max_it, double tol):
		world(world), A(world.ni*world.nj),
		max_solver_it(max_it), tolerance(tol) {
			buildMatrix();
			init();
		}

	/*computes electric field = -gradient(phi)*/
	void computeGradient();

	/*builds the "A" matrix for linear potential solver*/
    void buildMatrix();

    void init();

    //advance fields
    void advance();

    Field3 curl(Field3 &f, bool inner=false);
	Field div(Field3 &f);
	Field3 grad(Field &f) {return grad(f,world);}
	static Field3 grad(Field &f, World &world);

	static bool solveGSLinear(Matrix &A, dvector &x, dvector &b, int max_solver_it, double tolerance);

protected:
	World &world;
    Matrix A;				//system matrix for the linear equation


    enum NodeType {REG,NEUMANN,DIRICHLET};
    std::vector<NodeType> node_type;	//flag for different node types

	unsigned max_solver_it;	//maximum number of solver iterations
	double tolerance;		//solver tolerance

	/*linear GS solver for Ax=b system*/
	bool linearSolveGS(Matrix &A, dvector &x, dvector &b);
	void applyBoundaries();
};
#endif
