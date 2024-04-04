//ch1.cpp 
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>	

/*constants*/
namespace Const
{
	const double QE = 1.602176565e-19;	  // C, electron charge
	const double EPS_0 = 8.85418782e-12;  // C/V/m, vacuum perm.
	const double ME = 9.10938215e-31;	  // kg, electron mass
};

using namespace std;
using namespace Const;
using dvector = vector<double>;

/*function prototypes*/
bool outputCSV(double x0, double dx, const dvector &phi, const dvector &rho, const dvector &ef);
void solvePotentialDirect(double dx, dvector &phi, const dvector &rho);
bool solvePotentialGS(double dx, dvector &phi, const dvector &rho, int max_it=5000);
void computeEF(double dx, dvector &ef, const dvector &phi, bool second_order=true);
double gather(double li, dvector &field);

/*converts physical position x to a logical coordinate l*/
double XtoL(double x, double dx, double x0=0) {return (x-x0)/dx;}

/*main*/
int main()
{
	const int ni = 21;	//number of nodes
	const double x0 = 0;	//origin
	const double xd = 0.1;	//opposite end	
	double dx = (xd-x0)/(ni-1);	//node spacing

	dvector rho(ni,QE*1e12);
	dvector ef(ni);
	dvector phi(ni);
	
	//solve potential
	solvePotentialGS(dx,phi,rho);
	//solvePotentialDirect(dx,phi,rho);	//alternate solver

	//compute electric field
	computeEF(dx, ef, phi, true);

	//generate a test electron
	double m = ME;
	double q = -QE;
	double x = 4*dx;	//four cells from left edge
	double v = 0;		//stationary
	
	double dt = 1e-10;	//timestep

	//velocity rewind
	double li = XtoL(x,dx);
	double ef_p = gather(li,ef);
	v -= 0.5*(q/m)*ef_p*dt;
	
	//save initial potential for PE calculation
	double phi_max = phi[0];
	for (int i=1;i<ni;i++) 
		if (phi[i]>phi_max) phi_max = phi[i];
	
	//open file for particle trace
	ofstream out("trace.csv");
	if (!out) {cerr<<"Failed to open trace file"<<endl;return -1;}
	out<<"time,x,v,KE,PE\n";
	double x_old = x;
	
	//particle loop
	for (int ts=1;ts<=4000;ts++)
	{	
		//sample mesh data at particle position
		double li = XtoL(x,dx);
		double ef_p = gather(li,ef);
			
		//integrate velocity and position		
		x_old = x;
		v += (q/m)*ef_p*dt;
		x += v*dt;

		double phi_p = gather(XtoL(0.5*(x+x_old),dx),phi);
		double ke = 0.5*m*v*v/QE;			//KE in eV
		double pe = q*(phi_p-phi_max)/QE;	//PE in eV

		//write to file
		out<<ts*dt<<","<<x<<","<<v<<","<<ke<<","<<pe<<"\n";	
	
		if (ts==1 || ts%1000==0) //screen output every 1000 timesteps
			cout<<"ts: "<<ts<<", x:"<<x<<", v:"<<v<<", KE:"<<ke<<", PE:"<<pe<<"\n";		
	}
	
	//ouput to a CSV file for plotting
	outputCSV(ni,dx,phi,rho,ef);	
		
	return 0;	//normal exit
}

/*outputs the given fields to a CSV file, returns true if success*/
bool outputCSV(double x0, double dx, const dvector &phi, const dvector &rho, const dvector &ef)
{
	ofstream out("results.csv");	//open file for writing
	if (!out)
	{
		cerr<<"Could not open output file!"<<endl; 
		return false;
	}
	
	out<<"x,phi,rho,ef\n";		//write header
	for (size_t i=0;i<phi.size();i++)
	{
		out<<x0+i*dx; //write i-th position
		out<<","<<phi[i]<<","<< rho[i]<<","<<ef[i]; //write values
		out<<"\n";	//new line, not using endl to avoid buffer flush
	}
	
	return true; //file closed automatically here
}

/*solves Poisson's equation with Dirichlet boundaries using the Thomas algorithm*/
void solvePotentialDirect(double dx, dvector &phi, const dvector &rho)
{
	int ni = phi.size();	//number of mesh nodes		
	dvector a(ni);			//allocate memory for the matrix coefficients
	dvector b(ni);			
	dvector c(ni);
	dvector d(ni);	
			
	//set coefficients
	for (int i=0;i<ni;i++)
	{
		if (i==0 || i==ni-1)	//Dirichlet boundary
		{
			b[i] = 1;		//1 on the diagonal
			d[i] = 0;		//0 V 
		}
		else
		{
			a[i] = 1/(dx*dx);
			b[i] = -2/(dx*dx);
			c[i] = 1/(dx*dx);
			d[i] = -rho[i]/EPS_0; 
		}
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
	phi[ni-1] = d[ni-1];
	for (int i=ni-2;i>=0;i--)
	{
		phi[i] = d[i] - c[i]*phi[i+1];
	}
}

/* solves potential using the Gauss Seidel Method, returns true if converged*/
bool solvePotentialGS(double dx, dvector &phi, const dvector &rho, int max_it)
{
	double L2;
	double dx2 = dx*dx;		//precompute dx*dx
	const double w = 1.4;
	int ni = phi.size();	//number of mesh nodes
	
	/*solve potential*/
	for (int solver_it=0;solver_it<max_it;solver_it++)
	{
		phi[0] = 0;			//dirichlet boundary on left
		phi[ni-1] = 0;		//dirichlet boundary on right
		
		/*Gauss Seidel method, phi[i-1]-2*phi[i]+phi[i+1] = -dx^2*rho[i]/eps_0*/
		for (int i=1;i<ni-1;i++)
		{
			double g = 0.5*(phi[i-1] + phi[i+1] + dx2*rho[i]/EPS_0);
			phi[i] = phi[i] + w*(g-phi[i]);	//SOR
		}
			
		/*check for convergence*/
		if (solver_it%50==0)
		{
			double sum = 0;
			
			//internal nodes, automatically satisfied on Dirichlet boundaries
			for (int i=1;i<ni-1;i++)
			{
				double R = -rho[i]/EPS_0 - (phi[i-1] - 2*phi[i] + phi[i+1])/dx2;
				sum+=R*R;
			}
			L2 = sqrt(sum)/ni;
			if (L2<1e-6) 
			{
				cout<<"Gauss-Seidel converged after "<<solver_it<<" iterations"<<endl;
				return false;
			}
		}
	}
	cout<<"Gauss-Seidel failed to converge, L2="<<L2<<endl;
	return true;	
}

/* computes electric field by differentiating potential*/
void computeEF(double dx, dvector &ef, const dvector &phi, bool second_order)
{
	int ni = phi.size();	//number of mesh nodes

	//central difference on internal nodes
	for (int i=1;i<ni-1;i++)
		ef[i] = -(phi[i+1]-phi[i-1])/(2*dx);	

	//boundaries
	if (second_order)
	{
		ef[0] = (3*phi[0]-4*phi[1]+phi[2])/(2*dx);
		ef[ni-1] = (-phi[ni-3]+4*phi[ni-2]-3*phi[ni-1])/(2*dx);
	}
	else	//first order
	{
		ef[0] = (phi[0]-phi[1])/dx;
		ef[ni-1] = (phi[ni-2]-phi[ni-1])/dx;
	}
}

/*uses linear interpolation to evaluate f at li*/
double gather(double li, dvector &f)
{
	int i = (int)li;
	double di = li-i;
	return f[i]*(1-di) + f[i+1]*(di);
}

