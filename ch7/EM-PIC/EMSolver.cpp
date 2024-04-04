#include <math.h>
#include <iostream>
#include "World.h"
#include "EMSolver.h"
#include "Field.h"

using namespace std;
using namespace Const;

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
/*returns sum of v1[i]*v2[i]*/
double vec::dot(dvector v1, dvector v2)
{
	double dot = 0;
	size_t nu = v1.size();
	for (size_t j=0;j<nu;j++)
		dot+=v1[j]*v2[j];
	return dot;
}

/*returns l2 norm*/
double vec::norm(dvector v)
{
	double sum = 0;
	int nu = v.size();
	for (int j=0;j<nu;j++)
		sum+=v[j]*v[j];
	return sqrt(sum/nu);
}

/** converts 3D field to a 1D vector*/
dvector vec::deflate(Field &f3)
{
	dvector r(f3.ni*f3.nj);
	for (int i=0;i<f3.ni;i++)
		  for (int j=0;j<f3.nj;j++)
			 r[f3.U(i,j)] = f3[i][j];
	return r;
}

/** converts 1D vector to 3D field*/
void vec::inflate(dvector &d1, Field& f3)
{
	for (int i=0;i<f3.ni;i++)
		for (int j=0;j<f3.nj;j++)
			f3[i][j] = d1[f3.U(i,j)];
}

/** converts 1D vector to 3D field*/
Field vec::inflate(dvector &d1,int ni, int nj)
{
	Field f3(ni,nj);
	for (int i=0;i<ni;i++)
		for (int j=0;j<nj;j++)
			f3[i][j] = d1[f3.U(i,j)];
	return f3;
}

//constructs the coefficient matrix
void EMSolver::buildMatrix()
{
	double2 dh = world.getDh();
	double idx = 1.0/dh[0];
	double idy = 1.0/dh[1];
	double idx2 = idx*idx;	/*1/(dx*dx)*/
	double idy2 = idy*idy;
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
            if (i==0) {A(u,u)=idx;A(u,u+1)=-idx;}
            else if (i==ni-1) {A(u,u)=idx;A(u,u-1)=-idx;}
            else if (j==0) {A(u,u)=idy;A(u,u+ni)=-idy;}
            else if (j==nj-1) {A(u,u)=idy;A(u,u-ni)=-idy;}
            else {
            	//standard internal stencil
            	A(u,u-ni) = idy2;
            	A(u,u-1) = idx2;
            	A(u,u) = -2.0*(idx2+idy2);
            	A(u,u+1) = idx2;
            	A(u,u+ni) = idy2;
            	node_type[u] = REG;	//regular internal node
            }
        }
}


/*solves non-linear Poisson equation using Gauss-Seidel*/
bool EMSolver::linearSolveGS(Matrix &A, dvector &x, dvector &b)
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

//initial electric field
void EMSolver::init() {
	Field b3 = -(1/Const::EPS_0)*world.rho;

	//using solver operation on a 1D vector
	dvector b = vec::deflate(b3);
	dvector phi_s = vec::deflate(world.phi);
	//correct Dirichlet values
	for (int u=0;u<A.nu;u++)
		if (node_type[u]==NodeType::DIRICHLET) b[u] = phi_s[u];
	linearSolveGS(A,phi_s,b);
	world.phi = vec::inflate(phi_s, world.ni, world.nj);

	world.E = -1*grad(world.phi);
	applyBoundaries();

	// Faraday law
	Field3 curlE = curl(world.E);
	world.B += 0.5*world.getDt()*curlE;		//rewind by 0.5*dt

}

/*computes gradient(f) using 2nd order differencing*/
Field3 EMSolver::grad(Field &f, World &world)
{
	Field3 gf(f.ni,f.nj);	//this should be optimized out

	double2 dh = world.getDh();
	double dx = dh[0];
	double dy = dh[1];
	
	for (int i=0;i<world.ni;i++)
		for (int j=0;j<world.nj;j++)
		{
			/*x component*/
			if (i==0)
				gf[i][j][0] = -(-3*f[i][j]+4*f[i+1][j]-f[i+2][j])/(2*dx);	/*forward*/
			else if (i==world.ni-1)
				gf[i][j][0] = -(f[i-2][j]-4*f[i-1][j]+3*f[i][j])/(2*dx);	/*backward*/
			else
				gf[i][j][0] = -(f[i+1][j] - f[i-1][j])/(2*dx);	/*central*/

			/*y component*/
			if (j==0)
				gf[i][j][1] = -(-3*f[i][j] + 4*f[i][j+1]-f[i][j+2])/(2*dy);
			else if (j==world.nj-1)
				gf[i][j][1] = -(f[i][j-2] - 4*f[i][j-1] + 3*f[i][j])/(2*dy);
			else
				gf[i][j][1] = -(f[i][j+1] - f[i][j-1])/(2*dy);

			gf[i][j][2] = 0;	//no z component
		}
	return gf;
}

/* computes curl(f) for 2D problem with d()/dz=0 on staggered grid*/
Field3 EMSolver::curl(Field3 &f, bool inner) {
	int ni = f.ni;
	int nj = f.nj;

	if (inner) {ni++;nj++;}		//outer mesh is larger

	Field3 res(ni,nj);	//this should be optimized out
	double dx = world.getDh()[0];
	double dy = world.getDh()[1];

	int o = inner?1:0;

	for (int i=0;i<f.ni-1;i++)
		for (int j=0;j<f.nj-1;j++) {
			double3 fs = 0.5*(f[i][j]+f[i+1][j]);
			double3 fw = 0.5*(f[i+1][j]+f[i+1][j+1]);
			double3 fn = 0.5*(f[i][j+1]+f[i+1][j+1]);
			double3 fe = 0.5*(f[i][j]+f[i][j+1]);

			double dfz_dy = (fn[2]-fs[2])/dy;
			double dfy_dz = 0;
			double dfz_dx = (fw[2]-fe[2])/dx;
			double dfx_dz = 0;
			double dfy_dx = (fw[1]-fe[1])/dx;
			double dfx_dy = (fn[0]-fs[0])/dy;

			res[i+o][j+o][0] = dfz_dy - dfy_dz;
			res[i+o][j+o][1] = dfx_dz - dfz_dx;
			res[i+o][j+o][2] = dfy_dx - dfx_dy;
		}
	return res;
}

//evaluates divergence of f, only on internal nodes
Field EMSolver::div(Field3 &f) {
	Field div_f(f.ni, f.nj);

	double2 dh = world.getDh();
	double dx = dh[0];
	double dy = dh[1];

	for (int i=1;i<world.ni-1;i++)
		for (int j=1;j<world.nj-1;j++)
		{
			double dfx_dx = 0;
			double dfy_dy = 0;
			dfx_dx = (f[i+1][j][0] - f[i-1][j][0])/(2*dx);
			dfy_dy = (f[i][j+1][1] - f[i][j-1][1])/(2*dy);
			div_f[i][j] = dfx_dx + dfy_dy;	//no z component
		}
	return div_f;
}




void EMSolver::advance() {
	double dt = world.getDt();

	Field3 curlE = curl(world.E);
	world.B -= dt*curlE;		//B(k+0.5) = B(k-0.5) - curl[E(k)]*dt;

	// Ampere's law
	Field3 curlB = curl(world.B, true);
	world.E += dt*(Const::C*Const::C*curlB-(1/Const::EPS_0)*world.j);

	applyBoundaries();

}

void EMSolver::applyBoundaries() {

	//apply bc

    double d = 2;
    for (int i=0;i<world.ni;i++)
		for (int j=0;j<world.nj;j++)
		{

			double f=1,g=1;
			int L=world.ni-1;
			int M=world.nj-1;

			if (i<d)
				f = -i*i/(d*d)+2*i/d;
			else if (i>L-d) //if (x>(L-d))
				f = -i*i/(d*d) + 2*(L-d)*i/(d*d) - (L-d)*(L-d)/(d*d)+1;
			else f=1;

			if (j<d)
				g = -j*j/(d*d)+2*j/d;
			else if (j>M-d) //if (x>(L-d))
				g = -j*j/(d*d) + 2*(M-d)*j/(d*d) - (M-d)*(M-d)/(d*d)+1;
			else g= 1;
			world.E[i][j][0]*=f;
			world.E[i][j][1]*=g;

		}
}

/*solves non-linear Poisson equation using Gauss-Seidel*/
bool EMSolver::solveGSLinear(Matrix &A, dvector &x, dvector &b,
		int max_solver_it, double tolerance)
{
    double L2=0;			//norm
    bool converged= false;

    /*solve potential*/
    for (int it=0;it<max_solver_it;it++)
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


