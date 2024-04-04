/*
Steady-State Convection-Diffusion Equation solver
*/

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
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
	int ni = 81;

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
	
	//repeat for CDS
	double AcE = rho*u/(2*dx);
	double AcW = -rho*u/(2*dx);
	double AcP = -(AcW+AcE);

	//contribution from both terms
	double Aw = AdW + AcW;
	double Ap = AdP + AcP;
	double Ae = AdE + AcE;
		
	//set internal nodes
	for (int i=1;i<ni-1;i++) 
	{
		A.a[i] = Aw;
		A.b[i] = Ap;
		A.c[i] = Ae;
	}
	
	//initial value
	dvector phi(ni);
	phi[0] = b[0];
	phi[ni-1] = b[ni-1];

	//temporary vector for phi[k+1]
	dvector phi_new(ni);	
	
	//iterate using forward time
	double dt = 1e-2;
	
	//integrate solution in time
	for (int it=0;it<100;it++)
	{		
		//set only the non-boundary nodes
		for (int i=1;i<ni-1;i++)
		{
			//compute (A*phi) for node i
			double R = A.a[i]*phi[i-1] + 
				       A.b[i]*phi[i] + 
					   A.c[i]*phi[i+1];
			phi_new[i] = phi[i] - dt*R;
		}
		
		//copy down non-boundary nodes
		for (int i=1;i<ni-1;i++)
			phi[i] = phi_new[i];
		
		//open out file in the results folder
		stringstream ss;
		ss<<"results/ftcs_"<<setw(4)<<setfill('0')<<it<<".csv";		
		ofstream out(ss.str());
		out<<"x,phi"<<endl;
		for (int i=0;i<ni;i++)
			out<<L*i/(ni-1.0)<<","<<phi[i]<<"\n";
		out.close();		
	}
	cout<<"Done!"<<endl;
	return 0;
}		    
		
