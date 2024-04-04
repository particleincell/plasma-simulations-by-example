#include <math.h>
#include <iostream>
#include "World.h"
#include "PotentialSolver.h"
#include "Field.h"

using namespace std;
using namespace Const;

/*solves poisson equation with Boltzmann electrons using the Gauss-Seidel scheme*/

#include "PotentialSolver.h"
#include "Field.h"
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "World.h"

using namespace std;
using dvector = vector<double>;

//matrix-vector multiplication
dvector Matrix::operator*(dvector &v) {
	dvector r(nu);
	for (int u=0;u<nu;u++) {
		auto &row = rows[u];
		r[u] = 0;
		for (int i=0;i<nvals;i++){
			if (row.col[i]>=0) r[u]+=row.a[i]*v[row.col[i]];
			else break;	//end at the first -1
		}
	}
	return r;
}

//returns reference to A[r,c] element in the full matrix
double& Matrix::operator()(int r, int c){
	//find this entry
	auto &row = rows[r]; int v;
	for (v=0;v<nvals;v++)
	{
		if (row.col[v]==c) break;	//if found
		if (row.col[v]<0) {row.col[v]=c;   //set
						   break;}
	}
	assert(v!=nvals);	//check for overflow
	return row.a[v];
}

/*returns inverse of a diagonal preconditioner*/
Matrix Matrix::invDiagonal()
{
	Matrix M(nu);
	for (int r=0;r<nu;r++)	M(r,r) = 1.0/(*this)(r,r);
   return M;
}

/*subtracts diagonal matrix diag from A*/
Matrix Matrix::diagSubtract(dvector &P) {
	Matrix M(*this);	//make a copy
	for (int u=0;u<nu;u++) M(u,u)=(*this)(u,u)-P[u];
	return M;
}

//multiplies row r with vector x
double Matrix::multRow(int r, dvector &x){
	auto &row = rows[r];
	double sum=0;
	for (int i=0;i<nvals;i++)
	{
		if (row.col[i]>=0) sum+=row.a[i]*x[row.col[i]];
		else break;
	}
	return sum;
}


dvector operator-(const dvector &a, const dvector &b) {
	size_t nu = a.size();
	dvector r(nu);
	for (size_t u=0;u<nu;u++) r[u] = a[u]-b[u];
	return r;
}

dvector operator+(const dvector &a, const dvector &b) {
	size_t nu = a.size();
	dvector r(nu);
	for (size_t u=0;u<nu;u++) r[u] = a[u]+b[u];
	return r;
}

dvector operator*(const double s, const dvector &a) {
	size_t nu = a.size();
	dvector r(nu);
	for (size_t u=0;u<nu;u++) r[u] = s*a[u];
	return r;
}

/*vector math helper functions*/
namespace vec
{
	/*returns sum of v1[i]*v2[i]*/
	double dot(dvector v1, dvector v2)
	{
	    double dot = 0;
	    size_t nu = v1.size();
        for (size_t j=0;j<nu;j++)
            dot+=v1[j]*v2[j];
        return dot;
	}

	/*returns l2 norm*/
	double norm(dvector v)
	{
		double sum = 0;
		int nu = v.size();
        for (int j=0;j<nu;j++)
            sum+=v[j]*v[j];
		return sqrt(sum/nu);
	}

	/** converts 3D field to a 1D vector*/
	dvector deflate(Field &f3)
	{
		dvector r(f3.ni*f3.nj);
		for (int i=0;i<f3.ni;i++)
			  for (int j=0;j<f3.nj;j++)
				 r[f3.U(i,j)] = f3[i][j];
		return r;
	}

	/** converts 1D vector to 3D field*/
	void inflate(dvector &d1, Field& f3)
	{
		for (int i=0;i<f3.ni;i++)
			for (int j=0;j<f3.nj;j++)
				f3[i][j] = d1[f3.U(i,j)];
	}

};

//constructs the coefficient matrix
void PotentialSolver::buildMatrix()
{
	double2 dh = world.getDh();
	double idz = 1.0/dh[0];
	double idr = 1.0/dh[1];
	double idz2 = idz*idz;	/*1/(dz*dz)*/
	double idr2 = idr*idr;
	int ni = world.ni;
	int nj = world.nj;
	int nu = ni*nj;

	/*reserve space for node types*/
	node_type.reserve(nu);

	/*build matrix*/
    for (int j=0;j<nj;j++)
    	for (int i=0;i<ni;i++)
        {
            int u = world.U(i,j);
            A.clearRow(u);
            //dirichlet node?
			if (world.object_id[i][j]>0)
            {
                A(u,u)=1;	//set 1 on the diagonal
                node_type[u] = DIRICHLET;
                continue;
            }

			//Neumann boundaries
			node_type[u] = NEUMANN;		//set default
            if (i==0) {A(u,u)=idz;A(u,u+1)=-idz;}
            else if (i==ni-1) {A(u,u)=idz;A(u,u-1)=-idz;}
            else if (j==0) {A(u,u)=idr;A(u,u+ni)=-idr;}
            else if (j==nj-1) {A(u,u)=idr;A(u,u-ni)=-idr;}
            else {
            	//standard internal stencil
            	A(u,u-ni) = idr2 - 0.5*idr;  // 1/dr^2 - 1/(2*dr)
            	A(u,u-1) = idz2;
            	A(u,u) = -2.0*(idz2+idr2);
            	A(u,u+1) = idz2;
            	A(u,u+ni) = idr2 + 0.5*idr;	// 1/dr^2 + 1/(2*dr)
            	node_type[u] = REG;	//regular internal node
            }
        }
	solveQN();
}

/*quasi-neutral potential solver*/
bool PotentialSolver::solveQN()
{
	Field& phi = world.phi;
    Field& rhoi = world.rho;
    double rho0 = n0*QE;
    double rho_ratio_min = 1e-6;

    for (int i=0;i<world.ni;i++)
        for (int j=0;j<world.nj;j++)
        {
            if (world.object_id[i][j]>0) continue; /*skip Dirichlet nodes*/

			double rho_ratio = rhoi[i][j]/rho0;
			if (rho_ratio<rho_ratio_min) rho_ratio=rho_ratio_min;
			phi[i][j] = phi0 + Te0*log(rho_ratio);
        }
    return true;
}

/*Newton Raphson solver for a nonlinear system, using PCG for the linear solve	*/
bool PotentialSolver::solveNR()
{
	/*main NR iteration loop*/
	const int NR_MAX_IT=20;		/*maximum number of NR iterations*/
	const double NR_TOL = 1e-3;
	int nu = A.nu;

	Matrix J(nu);
	dvector P(nu);
	dvector y(nu);
	dvector x = vec::deflate(world.phi);
	dvector b = vec::deflate(world.rho);

	/*set RHS to zero on boundary nodes (zero electric field)
      and to existing potential on fixed nodes */
    for (int u=0;u<nu;u++)
    {
		if (node_type[u]==NEUMANN) b[u] = 0;			/*neumann boundary*/
        else if (node_type[u]==DIRICHLET) b[u] = x[u];	/*dirichlet boundary*/
        else b[u] = -b[u]/EPS_0;            /*regular node*/
    }

	double norm;
	bool converged=false;
	for(int it=0;it<NR_MAX_IT;it++)
	{
		/*compute F by first subtracting the linear term */
		dvector F = A*x-b;

		/*subtract b(x) on regular nodes*/
		for (int n=0;n<nu;n++)
			if (node_type[n]==REG)	/*regular nodes*/
				F[n] -= QE*n0*exp((x[n]-phi0)/Te0)/EPS_0;

		/*Compute P, diagonal of d(bx)/dphi*/
		for (int n=0;n<nu;n++)
		{
			if (node_type[n]==REG)
				P[n] = n0*QE/(EPS_0*Te0)*exp((x[n]-phi0)/Te0);
		}

		/*Compute J = A-diag(P)*/
		Matrix J = A.diagSubtract(P);

		/*solve Jy=F*/
        bool lin_converged = false;
		if (solver_type == PCG) 
			lin_converged = linearSolvePCG(J,y,F);
		if (!lin_converged)	//call GS if PCG failed to converge, or if not running with PCG
			linearSolveGS(J,y,F); 	   

		/*clear any numerical noise on Dirichlet nodes*/
		for (int u=0;u<nu;u++)
			if (node_type[u]==DIRICHLET) y[u]=0;

		/*x=x-y*/
		x = x-y;

		norm=vec::norm(y);
		//cout<<"NR norm: "<<norm<<endl;

		if (norm<NR_TOL)
		{
			converged=true;
			break;
		}
	}

	if (!converged)
		cout<<"NR+PCG failed to converge, norm = "<<norm<<endl;

	/*convert to 3d data*/
	vec::inflate(x,world.phi);
	return converged;
}

/*PCG solver for a linear system Ax=b*/
bool PotentialSolver::linearSolvePCG(Matrix &A, dvector &x, dvector &b)
{
	bool converged= false;

	double l2 = 0;
	Matrix M = A.invDiagonal(); //inverse of Jacobi preconditioner

	/*initialization*/
	dvector g = A*x-b;
	dvector s = M*g;
	dvector d = -1*s;

	for (unsigned it=0;it<max_solver_it;it++)
	{
		dvector z = A*d;
		double alpha = vec::dot(g,s);
		double beta = vec::dot(d,z);

		x = x+(alpha/beta)*d;
		g = g+(alpha/beta)*z;
		s = M*g;

		beta = alpha;
		alpha = vec::dot(g,s);

		d = (alpha/beta)*d-s;
		l2 = vec::norm(g);
		if (l2<tolerance) {converged=true;break;}
	}

	if (!converged)	cerr<<"PCG failed to converge, norm(g) = "<<l2<<endl;
    return converged;
}

/*solves non-linear Poisson equation using Gauss-Seidel*/
bool PotentialSolver::linearSolveGS(Matrix &A, dvector &x, dvector &b)
{
    double L2=0;			//norm
    bool converged= false;

    /*solve potential*/
    for (unsigned it=0;it<max_solver_it;it++)
    {
 		for (int u=0;u<A.nu;u++)
		{
			double S = A.multRow(u,x)-A(u,u)*x[u]; //multiplication of non-diagonal terms
 			double phi_new = (b[u]- S)/A(u,u);

 			/*SOR*/
            x[u] = x[u] + 1.*(phi_new-x[u]);
		}

		 /*check for convergence*/
		 if (it%25==0)
		 {
			 dvector R = A*x-b;
			 L2 = vec::norm(R);
			 if (L2<tolerance) {converged=true;break;}
		}
    }

    if (!converged) cerr<<"GS failed to converge, L2="<<L2<<endl;
    return converged;
}


/*computes electric field = -gradient(phi) using 2nd order differencing*/
void PotentialSolver::computeEF()
{
	//grab references to data
	Field &phi = world.phi;
	Field3 &ef = world.ef;

	double2 dh = world.getDh();
	double dx = dh[0];
	double dy = dh[1];
	
	for (int i=0;i<world.ni;i++)
		for (int j=0;j<world.nj;j++)
		{
			/*x component*/
			if (i==0)
				ef[i][j][0] = -(-3*phi[i][j]+4*phi[i+1][j]-phi[i+2][j])/(2*dx);	/*forward*/
			else if (i==world.ni-1)
				ef[i][j][0] = -(phi[i-2][j]-4*phi[i-1][j]+3*phi[i][j])/(2*dx);	/*backward*/
			else
				ef[i][j][0] = -(phi[i+1][j] - phi[i-1][j])/(2*dx);	/*central*/

			/*y component*/
			if (j==0)
				ef[i][j][1] = -(-3*phi[i][j] + 4*phi[i][j+1]-phi[i][j+2])/(2*dy);
			else if (j==world.nj-1)
				ef[i][j][1] = -(phi[i][j-2] - 4*phi[i][j-1] + 3*phi[i][j])/(2*dy);
			else
				ef[i][j][1] = -(phi[i][j+1] - phi[i][j-1])/(2*dy);

			ef[i][j][2] = 0;	//no z component
		}
}
