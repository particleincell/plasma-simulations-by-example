/*
Example Vlasov code for 2 stream instability
see https://www.particleincell.com/2018/into-to-vlasov-solvers/
Written by Lubos Brieda
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <cstring>
#include <cmath>
#include <vector>
#include "Field.h"

using namespace std;

//storage of domain parameters
struct World
{
	double L;	//length of x
	double v_max;	//extends in the velocity space
	double dx,dv;	//cell spacing
	int ni,nj;	//number of nodes in x,y
	bool periodic = true;	//controls whether the world is periodic in x
	
	void setLimits(double L, double v_max) {this->L=L;this->v_max=v_max;}
	void setNodes(int N, int M) {ni=N; nj=2*M-1;dx=L/(ni-1);dv=2*v_max/(nj-1);} 
	double getX(int i) {return 0+i*dx;}
	double getV(int j) {return -v_max+j*dv;}
	
	//linear interpolation: a higher order scheme needed!
	double interp(Field &f, double x, double v)
	{
		double fi = (x-0)/dx;
		double fj = (v-(-v_max))/dv;		
		
		//periodic boundaries in i
		if (periodic)
		{
			if (fi<0) fi+=ni-1;		//-0.5 becomes ni-1.5, which is valid since i=0=ni-1
			if (fi>ni-1) fi-=ni-1;
		}
		else if (fi<0 || fi>=ni-1) return 0;
		
		//return zero if velocity less or more than limits
		if (fj<0 || fj>=nj-1) return 0;
		
		int i = (int)fi;
		int j = (int)fj;
		double di = fi-i;
		double dj = fj-j;
		
		double val = (1-di)*(1-dj)*f[i][j];
		if (i<ni-1) val+=(di)*(1-dj)*f[i+1][j];
		if (j<nj-1) val+=(1-di)*(dj)*f[i][j+1];
		if (i<ni-1 && j<nj-1) val+=(di)*(dj)*f[i+1][j+1];
		return val;		
	}
	
	/*makes values on left and right edge identical on periodic systems*/
	void applyBC(Field &f)
	{
		if (!periodic)	return;
		for (int j=0;j<nj;j++)
		{
			f[0][j] = 0.5*(f[0][j]+f[ni-1][j]);
			f[ni-1][j] = f[0][j];
		}
	}
	
};

//filter to eliminate garbage values like 1e-120
double filter(double a) {if (std::abs(a)<1e-20) return 0; else return a;}

/*saves the provided scalars and vectors to a file*/
void saveVTK(int time_step, World &world, map<string,Field*> scalars2D, map<string,dvector*> scalars1D)
{
	//generate file name
	stringstream ss;
	ss<<"results/vlasov";
	if (time_step>=0)
		ss<<"_"<<setw(6)<<setfill('0')<<time_step;
	ss<<".vti";	
	ofstream out(ss.str());

	out<<setprecision(4);

	out<<"<VTKFile type=\"ImageData\">\n";
	out<<"<ImageData Origin=\""<<0<<" "<<-world.v_max<<" "<<0<<"\"";
	out<<" Spacing=\""<<world.dx<<" "<<world.dv<<" "<<1<<"\"";
	out<<" WholeExtent=\""<<0<<" "<<world.ni-1<<" "<<0<<" "<<world.nj-1<<" "<<0<<" "<<0<<"\">\n";
	out<<"<Piece Extent=\""<<0<<" "<<world.ni-1<<" "<<0<<" "<<world.nj-1<<" "<<0<<" "<<0<<"\">\n";
	out<<"<PointData>\n";
		
	//user vars, p.first is the string key, p.second is the double* pointer to data
	for (pair<string,Field*> p : scalars2D)
	{
		//p.first is the string key, p.second is the double* pointer to data
		out<<"<DataArray Name=\""<<p.first<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		Field &f = *p.second;
		
		for (int j=0;j<world.nj;j++)
			for (int i=0;i<world.ni;i++)
				out<<setprecision(4)<<filter(f[i][j])<<" ";					
			out<<"\n";
		out<<"</DataArray>\n";
	}

	for (pair<string,dvector*> p : scalars1D)
	{
		//p.first is the string key, p.second is the double* pointer to data
		out<<"<DataArray Name=\""<<p.first<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		dvector &f = *p.second;
		
		for (int j=0;j<world.nj;j++)
			for (int i=0;i<world.ni;i++)
				out<<setprecision(4)<<filter(f[i])<<" ";					
			out<<"\n";
		out<<"</DataArray>\n";
	}

	out<<"</PointData>\n";
	out<<"</Piece>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";

	out.close();
}

/*solves Poisson's equation with Dirichlet boundaries using the direct Thomas algorithm
and returns the electric field*/
void solvePoissonsEquationGS(World &world, dvector &b, dvector &phi, dvector &E)
{
	double dx2 = world.dx*world.dx;
	int ni=world.ni;
	double norm;
	double const tol = 1e-3;

	for (int i=0;i<ni;i++) phi[i]=0;		//clear data

	for (int it=0;it<10000;it++)
	{
		/*periodic boundary conditions*/
		phi[0] = 0.5*(phi[ni-2]+phi[1]-dx2*b[0]);
		
		for (int i=1;i<ni-1;i++)
		{
			double g = 0.5*(phi[i-1]+phi[i+1]-dx2*b[i]);
			phi[i] = phi[i] + 1.4*(g-phi[i]);
		}

		phi[ni-1] = 0.5*(phi[ni-2]+phi[1]-dx2*b[ni-1]);
			
		/*check for convergence*/
		if (it%50 == 0)
		{
			double R_sum = 0;
			for (int i=1;i<ni-1;i++)
			{
				double dR =  (phi[i-1]-2.*phi[i]+phi[i+1])/dx2 - b[i];
				R_sum += dR*dR;				
			}
			//periodic boundaries
			double dR = (phi[ni-2]-2*phi[0]+phi[1])/dx2 - b[0];
			R_sum += dR*dR;
			dR = (phi[ni-2]-2*phi[ni-1]+phi[1])/dx2 - b[ni-1];
			R_sum += dR*dR;
			
			norm = sqrt(R_sum/ni);
			//cout<<norm<<endl;
			if (norm<tol)
				break;
		}
	}
	
	if (norm>tol)
		cout<<"GS failed to converge, norm = "<<norm<<endl;
		
	/*set periodic boundary*/	
	phi[0] = 0.5*(phi[0]+phi[ni-1]);
    phi[ni-1] = phi[0];
	
	//compute electric field
	for (int i=1;i<ni-1;i++) 
		E[i] = -(phi[i+1]-phi[i-1])/(2*world.dx);
	E[0] = E[ni-1] = -(phi[1]-phi[ni-2])/(2*world.dx);
}

int main()
{
	//constants and parameters
	const double pi = acos(-1.0);		//pi

	//create a variable of type World
	World world;	
	world.setLimits(10,5);
	world.setNodes(101,51);
	int ni = world.ni;
	int nj = world.nj;
	double dx = world.dx;
	double dv = world.dv;
	double dt = 1/8.0;
	
	cout<<"dx: "<<dx<<" "<<"dv: "<<dv<<endl;
			
	Field f(ni,nj); 	//f
	Field fs(ni,nj);	//fs
	Field fss(ni,nj);	//fss
	dvector ne(ni);		//number density
	dvector b(ni);		//Poisson solver RHS, -rho=(ne-1) since ni=1 is assumed
	dvector E(ni);		//electric field
	dvector phi(ni);	//potential
	
	//map is a list of keys and corresponding values for output
	map<string,Field*> scalars2D; 
	map<string,dvector*> scalars1D; 
	
	scalars2D["f"] = &f;	
	scalars1D["ne"] = &ne;
	scalars1D["E"] = &E;
	
	world.periodic = true;

	//set initial distribution
	for (int i=0;i<ni;i++)
		for (int j=0;j<nj;j++)
		{
			double x = world.getX(i);
			double v = world.getV(j);			
						
			double vth2 = 0.02;
			double vs1 = 1.6;
			double vs2 = -1.4;
			f[i][j] = 0.5/sqrt(vth2*pi)*exp(-(v-vs1)*(v-vs1)/vth2);
			f[i][j] += 0.5/sqrt(vth2*pi)*exp(-(v-vs2)*(v-vs2)/vth2)*(1+0.02*cos(3*pi*x/world.L));
		}
	
	//set some constant e field
	for (int i=0;i<ni;i++)
		E[i] = 0;
	
	int it;
	
	//main loop
	for (it=0;it<=1000;it++)
	{
		if (it%100==0)	cout<<it<<endl;
		if (it%5==0) saveVTK(it,world,scalars2D,scalars1D);

		//compute f*
		for (int i=0;i<ni;i++)
			for (int j=0;j<nj;j++)
			{
				double v = world.getV(j);
				double x = world.getX(i);
				
				fs[i][j] = world.interp(f,x-v*0.5*dt,v);								
			}
		
		world.applyBC(fs);
				
		//compute number density by integrating f with the trapezoidal rule		
		for (int i=0;i<ni;i++)
		{
			ne[i] = 0;
			for (int j=0;j<nj-1;j++)
				ne[i]+=0.5*(fs[i][j+1]+fs[i][j])*dv;
		}
		
		//compute the right hand side, -rho = (ne-1)
		for (int i=0;i<ni;i++)	b[i] = ne[i]-1;		
		b[0] = 0.5*(b[0]+b[ni-1]);
		b[ni-1] = b[0];

		//solution of the Poisson's equation
		solvePoissonsEquationGS(world,b,phi,E);
		
		//compute f**
		for (int i=0;i<ni;i++)
			for(int j=0;j<nj;j++)
			{
				double v = world.getV(j);
				double x = world.getX(i);
				fss[i][j] = world.interp(fs,x,v+E[i]*dt);				
			}
		
		world.applyBC(fss);
		
		//compute f(n+1)
		for (int i=0;i<ni;i++)
			for(int j=0;j<nj;j++)
			{
				double v = world.getV(j);
				double x = world.getX(i);
				f[i][j] = world.interp(fss,x-v*0.5*dt,v);
			}		
		
		world.applyBC(f);
	}
				
	saveVTK(it,world,scalars2D,scalars1D);
		
	return 0;
}
