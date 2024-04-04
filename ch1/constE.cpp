#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

/*constants*/
#define QE 1.602176565e-19		// C, electron charge
#define EPS_0 8.85418782e-12  	// F/m, vacuum permittivity
#define ME 9.10938215e-31		// kg, electron mass

/*main*/
int main()
{
	//load particle
	double m = ME;	//particle mass
	double q = -QE;	//particle charge
	double x = 0;	//initial position
	double v = 0;	//stationary
	
	//simulation parameters
	double dt = 1e-9;	//timestep
	double E = -100;	//electric field

	//open particle trace file
	ofstream out("trace.csv");
	if (!out) {cerr<<"Failed to open file"<<endl;return -1;}
	out<<"t,x,v,x2,v2,x3,v3\n";
		
	//particle loop
	for (int it=0;it<20;it++)
	{
		//write trace data
		out<<it*dt<<","<<x<<","<<v<<"\n";
		
		//compute future position and velocity
		x = x + v*dt; 
		v = v + (q/m)*E*dt;
	}
	
	return 0;	//normal exit
}
