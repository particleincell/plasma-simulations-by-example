/*
Steady-State Convection-Diffusion Equation solver
*/

#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>

using namespace std;
using dvector = vector<double>;

// simple tri diagonal matrix
struct TriD {
	TriD(int ni) {
		a.reserve(ni);
		b.reserve(ni);
		c.reserve(ni);
	}

	dvector a,b,c;	//coefficients for left, main, and right diagonals
};

/*solves Ax=b using the Thomas algorithm*/
dvector solveDirect(TriD &A, dvector &rhs)
{
	int ni = rhs.size();	//number of mesh nodes		
	dvector a(ni);			//reserve memory space
	dvector b(ni);			
	dvector c(ni);
	dvector d(ni);	
			
	//copy coefficients
	for (int i=0;i<ni;i++)
	{
		a[i] = A.a[i];
		b[i] = A.b[i];
		c[i] = A.c[i];
		d[i] = rhs[i];
	}
	
	//initialize
	c[0] = c[0]/b[0];
	d[0] = d[0]/b[0];
	
	//forward step
	for (int i=1;i<ni;i++)
	{
		if (i<ni-1)
			c[i] = c[i]/(b[i]-a[i]*c[i-1]);
		
		d[i] = (d[i] - a[i]*d[i-1])/(b[i] - a[i]*c[i-1]);
	}
	
	//backward substitution
	dvector x(ni);
	x[ni-1] = d[ni-1];
	for (int i=ni-2;i>=0;i--)
	{
		x[i] = d[i] - c[i]*x[i+1];
	}
	return x;
}

int main() 
{
	//simulation inputs
	double rho = 1;
	double u = 1;
	double D = 0.02;
	double phi0 = 0;
	double phiL = 1;

	//set domain parameters
 	double L = 1;
	int ni = 41;

	//write out analytical solution
	ofstream out_th("theory.csv");
	out_th<<"x,phi_th"<<endl;
	for (int i=0;i<501;i++) {
		double x = L*i/500.0;
		double Pe = rho*u*L/D;
		double phi_true = phi0 + ((exp(x*Pe/L)-1)/ (exp(Pe)-1)*(phiL-phi0));
		out_th<<x<<","<<phi_true<<"\n";
	}

	//set matrix
	TriD A(ni);
	dvector b(ni);
	
	//dirichlet condition on left and right
	A.b[0] = 1;
	b[0] = phi0;
	A.b[ni-1] = 1;
	b[ni-1] = phiL;

	//assuming uniform spacing
	double dx = L/(ni-1);

	//diffusive term
	double AdW = -D/(dx*dx);
	double AdE = -D/(dx*dx);
	double AdP = -(AdE + AdW);

	//upwind scheme for the convective derivative
	double AcE = min(rho*u,0.)/dx;
	double AcW = -max(rho*u,0.)/dx;
	double AcP = -(AcE+AcW);

	//contribution from both terms
	double Aw = AdW + AcW;
	double Ap = AdP + AcP;
	double Ae = AdE + AcE;

	//set internal nodes
	for (int i=1;i<ni-1;i++)
	{
		A.a[i] = Aw;  //A[i,i-1]
		A.b[i] = Ap;  //A[i,i]
		A.c[i] = Ae;  //A[i,i+1]
	}

	//obtain the solution    
	dvector phi_uds = solveDirect(A, b);
	
	//repeat for CDS
	AcE = rho*u/(2*dx);
	AcW = -rho*u/(2*dx);
	AcP = -(AcW+AcE);

	//contribution from both terms
	Aw = AdW + AcW;
	Ap = AdP + AcP;
	Ae = AdE + AcE;
		
	//set internal nodes
	for (int i=1;i<ni-1;i++) 
	{
		A.a[i] = Aw;
		A.b[i] = Ap;
		A.c[i] = Ae;
	}

	//obtain the solution    
	dvector phi_cds = solveDirect(A,b);

	//output results
	ofstream out("numerical41.csv");
	out<<"x,phi_uds,phi_cds"<<endl;
	for (int i=0;i<ni;i++) {
		double x = L*i/(ni-1.0);
		out<<x<<","<<phi_uds[i]<<","<<phi_cds[i]<<"\n";
	}
	return 0;
}		    
		
