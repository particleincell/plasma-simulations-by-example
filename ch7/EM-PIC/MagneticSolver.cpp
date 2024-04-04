#include <math.h>
#include <iostream>
#include "World.h"
#include "MagneticSolver.h"
#include "Field.h"

using namespace std;
using namespace Const;



//constructs the coefficient matrix
void MagneticSolver::buildMatrix()
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
	phi_m.reserve(nu);

	//get sphere center
	double3 orig_lc = world.XtoL(world.getCircleOrig());
	sphere_orig_u = world.U((int)orig_lc[0], (int)orig_lc[1]);

	/*solve potential*/
	for (int j=0;j<nj;j++)
		for (int i=0;i<ni;i++)
		{
			int u = world.U(i,j);
			A.clearRow(u);

			//if sphere center
			if (u==sphere_orig_u)
			{
				node_type[u] = DIRICHLET;
				A(u,u)=1;
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

dvector MagneticSolver::computeDivergence(Field3 &M)
{
	double3 dh = world.getDh();
	double dx = dh[0];
	double dy = dh[1];
	dvector div(A.nu);

	for (int i=0;i<world.ni;i++)
		for (int j=0;j<world.nj;j++)
		{
			double divX = 0, divY=0;

			/*x component*/
			if (i==0)
				divX = (-3*M[i][j][0]+4*M[i+1][j][0]-M[i+1][j][0])/(2*dx);	/*forward*/
			else if (i==world.ni-1)
				divX = (M[i-2][j][0]-4*M[i-1][j][0]+3*M[i][j][0])/(2*dx);	/*backward*/
			else
				divX = (M[i+1][j][0] - M[i-1][j][0])/(2*dx);	/*central*/

			/*y component*/
			if (j==0)
				divY = (-3*M[i][j][1] + 4*M[i][j+1][1]-M[i][j+2][1])/(2*dy);
			else if (j==world.nj-1)
				divY = (M[i][j-2][1] - 4*M[i][j-1][1] + 3*M[i][j][1])/(2*dy);
			else
				divY = (M[i][j+1][1] - M[i][j-1][1])/(2*dy);

			int u = j*world.ni+i;
			div[u] = divX+divY;
		}

	return div;
}

bool MagneticSolver::solve() {
	dvector b_m = computeDivergence(world.M);	//-rho_m=div(M)
	b_m[sphere_orig_u] = 0;		//Dirichlet solution at sphere center
	vec::inflate(b_m,world.b_m);

	bool conv = EMSolver::solveGSLinear(A, phi_m, b_m, 200000, 1e-6);
	vec::inflate(phi_m,world.phi_m); //set 3D solution for visualization

	world.H = -1*EMSolver::grad(world.phi_m, world);

	world.B = Const::MU_0*(world.H+world.M);


	return conv;
}


