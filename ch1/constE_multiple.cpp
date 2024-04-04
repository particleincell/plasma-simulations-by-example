/* integrates particle position using forward, backward, and central methods*/
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

/*constants*/
namespace Const {
	double QE = 1.602176565e-19;		// C, electron charge
	double EPS_0 = 8.85418782e-12;  	// F/m, vacuum permittivity
	double ME = 9.10938215e-31;		// kg, electron mass
}

using namespace Const;

/*main*/
int main()
{
	//load particle
	double m = ME;	//particle mass
	double q = -QE;	//particle charge
	double x = 0;	//initial position
	double v = 0;	//stationary
	
	double x2 = 0;	//initial position
	double v2 = 0;	//stationary
	
	double x3 = 0;
	double v3 = 0;
	
	//simulation parameters
	double dt = 1e-9;	//timestep
	double E = -100;	//electric field

	//open particle trace file
	ofstream out("trace.csv");
	if (!out) {cerr<<"Failed to open file"<<endl;return -1;}
	out<<"t,x,v,x2,v2,x3,v3\n";
		
	//particle loop
	for (int it=0;it<10;it++)
	{
		//write trace data
		out<<it*dt<<","<<x<<","<<v
				  <<","<<x2<<","<<v2
				  <<","<<x3<<","<<v3
				  <<"\n";
		
		//forward
		x = x + v*dt; 
		v = v + (q/m)*E*dt;
		
		//backward
		v2 = v2 + (q/m)*E*dt;
		x2 = x2 + v2*dt; 
		
		//central
		double v_new = v3 + (q/m)*E*dt;
		x3 = x3 + 0.5*(v_new+v3)*dt; 
		v3=v_new;
	}
		
	/*run again with smaller dt*/
	x = 0;	//initial position
	v = 0;	//stationary
	
	x2 = 0;	//initial position
	v2 = 0;	//stationary
	//open particle trace file
	ofstream out2("trace2.csv");
	if (!out2) {cerr<<"Failed to open file"<<endl;return -1;}
	out2<<"t,x,v,x2,v2\n";
    dt /= 2;
	//particle loop
	for (int it=0;it<10*2;it++)
	{
		//write trace data
		out2<<it*dt<<","<<x<<","<<v
				  <<","<<x2<<","<<v2
				  <<"\n";
		
		//compute future position and velocity
		x = x + v*dt; 
		v = v + (q/m)*E*dt;
		
		//compute future position and velocity
		v2 = v2 + (q/m)*E*dt;
		x2 = x2 + v2*dt; 
		
	}
	return 0;	//normal exit
}
