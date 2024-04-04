#include <math.h>
#include <iostream>
#include "World.h"
#include "MagneticSolver.h"
#include "PotentialSolver.h"
#include "Field.h"

using namespace std;
using namespace Const;

//constructs the coefficient matrix
void MagneticSolver::buildMatrix()
{
	double3 dh = world.getDh();
	double idx = 1.0/dh[0];
	double idy = 1.0/dh[1];
	double idz = 1.0/dh[2];
    double idx2 = idx*idx;	/*1/(dx*dx)*/
	double idy2 = idy*idy;
	double idz2 = idz*idz;
	int ni = world.ni;
	int nj = world.nj;
	int nk = world.nk;
	int nu = ni*nj*nk;

	/*reserve space for node types*/
	node_type.reserve(nu);
	phi_m.reserve(nu);

	//get sphere center
	double3 orig_lc = world.XtoL(world.getSphereOrig());
	sphere_orig_u = world.U((int)orig_lc[0], (int)orig_lc[1], (int)orig_lc[2]);

	/*solve potential*/
	for (int k=0;k<nk;k++)
        for (int j=0;j<nj;j++)
        	for (int i=0;i<ni;i++)
            {
                int u = world.U(i,j,k);
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
                else if (k==0) {A(u,u)=idz;A(u,u+ni*nj)=-idz;}
				else if (k==nk-1) {
					A(u,u)=idz;
					A(u,u-ni*nj)=-idz;}
                else {
                	//standard internal stencil
                	A(u,u-ni*nj) = idz2;
                	A(u,u-ni) = idy2;
                	A(u,u-1) = idx2;
                	A(u,u) = -2.0*(idx2+idy2+idz2);
                	A(u,u+1) = idx2;
                	A(u,u+ni) = idy2;
                	A(u,u+ni*nj) = idz2;
                	node_type[u] = REG;	//regular internal node
                }
            }
}

dvector MagneticSolver::computeDivergence(Field3 &M)
{
	double3 dh = world.getDh();
	double dx = dh[0];
	double dy = dh[1];
	double dz = dh[2];
	dvector div(A.nu);

	for (int i=0;i<world.ni;i++)
		for (int j=0;j<world.nj;j++)
			for (int k=0;k<world.nk;k++)
			{
				double divX = 0, divY=0, divZ=0;

				/*x component*/
				if (i==0)
					divX = (-3*M[i][j][k][0]+4*M[i+1][j][k][0]-M[i+1][j][k][0])/(2*dx);	/*forward*/
				else if (i==world.ni-1)
					divX = (M[i-2][j][k][0]-4*M[i-1][j][k][0]+3*M[i][j][k][0])/(2*dx);	/*backward*/
				else
					divX = (M[i+1][j][k][0] - M[i-1][j][k][0])/(2*dx);	/*central*/

				/*y component*/
				if (j==0)
					divY = (-3*M[i][j][k][1] + 4*M[i][j+1][k][1]-M[i][j+2][k][1])/(2*dy);
				else if (j==world.nj-1)
					divY = (M[i][j-2][k][1] - 4*M[i][j-1][k][1] + 3*M[i][j][k][1])/(2*dy);
				else
					divY = (M[i][j+1][k][1] - M[i][j-1][k][1])/(2*dy);

				/*z component*/
				if (k==0)
					divZ = (-3*M[i][j][k][2] + 4*M[i][j][k+1][2]-M[i][j][k+2][2])/(2*dz);
				else if (k==world.nk-1)
					divZ = (M[i][j][k-2][2] - 4*M[i][j][k-1][2]+3*M[i][j][k][2])/(2*dz);
				else
					divZ = (M[i][j][k+1][2] - M[i][j][k-1][2])/(2*dz);

				int u = k*world.ni*world.nj+j*world.ni+i;
				div[u] = divX+divY+divZ;
			}

	return div;
}

bool MagneticSolver::solve() {
	dvector b_m = computeDivergence(world.M);	//-rho_m=div(M)
	b_m[sphere_orig_u] = 0;		//Dirichlet solution at sphere center
	bool conv = PotentialSolver::solveGSLinear(A, phi_m, b_m, max_solver_it, tolerance);

	PotentialSolver::computeNegGradient(world.H,phi_m,world);
	world.B = Const::MU_0*(world.H+world.M);

	//visualization support
	vec::inflate(phi_m,world.phi_m); //set 3D solution for visualization
    analyticalSol();	//also compute analytical solution

	return conv;
}

void MagneticSolver::analyticalSol()
{
	Field &sol = world.phi_m_theory;
	double3 axis{0,0,1};		//z axis
	double a = world.getSphereRad();
	//get M0 by evaluating M at sphere center
	double3 sphere_x0 = world.getSphereOrig();
	double M0 = mag(world.M.gather(world.XtoL(sphere_x0)));

	for (int i=0;i<world.ni;i++)
		for (int j=0;j<world.nj;j++)
			for (int k=0;k<world.nk;k++)
			{
				double3 pos = world.pos(i,j,k);
				double3 ray = pos - sphere_x0;	//ray to the point from sphere origing
				double r = mag(ray);	//radius
				if (r>0) {
					double cos_theta = dot(ray,axis)/r;
					double rl = min(r,a);
					double rm = max(r,a);
					sol[i][j][k] = (1/3.0)*M0*a*a*rl/(rm*rm)*cos_theta;
				}
				else sol[i][j][k] = 0;
			}
}
