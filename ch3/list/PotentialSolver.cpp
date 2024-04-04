#include <math.h>
#include <iostream>
#include "World.h"
#include "PotentialSolver.h"
#include "Field.h"

using namespace std;
using namespace Const;

/*solves non-linear Poisson equation using Gauss-Seidel*/
bool PotentialSolver::solve()
{
    //references to avoid having to write world.phi
	Field &phi = world.phi;
    Field &rho = world.rho;		//rho contains only ion contribution

	//precompute 1/(dx^2)
    double3 dh = world.getDh();
    double idx2 = 1.0/(dh[0]*dh[0]);
    double idy2 = 1.0/(dh[1]*dh[1]);
    double idz2 = 1.0/(dh[2]*dh[2]);

    double L2=0;			//norm
    bool converged= false;

    /*solve potential*/
    for (unsigned it=0;it<max_solver_it;it++)
    {
 		for (int i=0;i<world.ni;i++)
        	for (int j=0;j<world.nj;j++)
            	for (int k=0;k<world.nk;k++)
                {
                    /*skip over solid (fixed) nodes = Dirichlet boundaries*/
                    if (world.object_id[i][j][k]>0) continue;

                    if (i==0)
                    	phi[i][j][k] = phi[i+1][j][k];
                    else if (i==world.ni-1)
                    	phi[i][j][k] = phi[i-1][j][k];
                    else if (j==0)
                    	phi[i][j][k] = phi[i][j+1][k];
                    else if (j==world.nj-1)
                    	phi[i][j][k] = phi[i][j-1][k];
                    else if (k==0)
                    	phi[i][j][k] = phi[i][j][k+1];
                    else if (k==world.nk-1)
                    	phi[i][j][k] = phi[i][j][k-1];
                    else {	//standard internal open node

                    	//evaluate electron density from the Boltzmann relationshp
                    	double ne = n0 * exp((phi[i][j][k]-phi0)/Te0);

                    	double phi_new = ((rho[i][j][k]-Const::QE*ne)/Const::EPS_0 +
                                        idx2*(phi[i-1][j][k] + phi[i+1][j][k]) +
                                        idy2*(phi[i][j-1][k]+phi[i][j+1][k]) +
                                        idz2*(phi[i][j][k-1]+phi[i][j][k+1]))/(2*idx2+2*idy2+2*idz2);

                    	/*SOR*/
                    	phi[i][j][k] = phi[i][j][k] + 1.4*(phi_new-phi[i][j][k]);
                    }
                }

		 /*check for convergence*/
		 if (it%25==0)
		 {
			double sum = 0;
			for (int i=0;i<world.ni;i++)
				for (int j=0;j<world.nj;j++)
					for (int k=0;k<world.nk;k++)
					{
						/*skip over solid (fixed) nodes*/
						if (world.object_id[i][j][k]>0) continue;

						double R = 0;
                    	if (i==0)
                    		R = phi[i][j][k] - phi[i+1][j][k];
                    	else if (i==world.ni-1)
                    		R = phi[i][j][k] - phi[i-1][j][k];
                    	else if (j==0)
                    		R = phi[i][j][k] - phi[i][j+1][k];
                    	else if (j==world.nj-1)
                    		R = phi[i][j][k] - phi[i][j-1][k];
                    	else if (k==0)
                    		R = phi[i][j][k] - phi[i][j][k+1];
                    	else if (k==world.nk-1)
                    		R = phi[i][j][k] - phi[i][j][k-1];
                    	else {
                    			//evaluate electron density from the Boltzmann relationshp
                    		    double ne = n0 * exp((phi[i][j][k]-phi0)/Te0);
                    			R = -phi[i][j][k]*(2*idx2+2*idy2+2*idz2) +
									(rho[i][j][k]-Const::QE*ne)/Const::EPS_0 +
									idx2*(phi[i-1][j][k] + phi[i+1][j][k]) +
									idy2*(phi[i][j-1][k]+phi[i][j+1][k]) +
									idz2*(phi[i][j][k-1]+phi[i][j][k+1]);
						}

						sum += R*R;
					}

			L2 = sqrt(sum/(world.ni*world.nj*world.nk));
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

	double3 dh = world.getDh();
	double dx = dh[0];
	double dy = dh[1];
	double dz = dh[2];

	for (int i=0;i<world.ni;i++)
		for (int j=0;j<world.nj;j++)
			for (int k=0;k<world.nk;k++)
			{
				/*x component*/
				if (i==0)
					ef[i][j][k][0] = -(-3*phi[i][j][k]+4*phi[i+1][j][k]-phi[i+2][j][k])/(2*dx);	/*forward*/
				else if (i==world.ni-1)
					ef[i][j][k][0] = -(phi[i-2][j][k]-4*phi[i-1][j][k]+3*phi[i][j][k])/(2*dx);	/*backward*/
				else
					ef[i][j][k][0] = -(phi[i+1][j][k] - phi[i-1][j][k])/(2*dx);	/*central*/

				/*y component*/
				if (j==0)
					ef[i][j][k][1] = -(-3*phi[i][j][k] + 4*phi[i][j+1][k]-phi[i][j+2][k])/(2*dy);
				else if (j==world.nj-1)
					ef[i][j][k][1] = -(phi[i][j-2][k] - 4*phi[i][j-1][k] + 3*phi[i][j][k])/(2*dy);
				else
					ef[i][j][k][1] = -(phi[i][j+1][k] - phi[i][j-1][k])/(2*dy);

				/*z component*/
				if (k==0)
					ef[i][j][k][2] = -(-3*phi[i][j][k] + 4*phi[i][j][k+1]-phi[i][j][k+2])/(2*dz);
				else if (k==world.nk-1)
					ef[i][j][k][2] = -(phi[i][j][k-2] - 4*phi[i][j][k-1]+3*phi[i][j][k])/(2*dz);
				else
					ef[i][j][k][2] = -(phi[i][j][k+1] - phi[i][j][k-1])/(2*dz);
			}
}
