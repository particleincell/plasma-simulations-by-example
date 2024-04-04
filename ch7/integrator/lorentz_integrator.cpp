//demo of particle integration with magnetic field with the Boris method
#include <fstream>
#include <iostream>
#include "Vec.h"

using namespace std;

//simple particle container
struct Particle {
	double3 x,v;	//position and velocity
	double q;		//charge
	double m;		//mass
	Particle(const double3 &x, const double3 &v, double q, double m):
	x{x}, v{v}, q{q}, m{m} {}
};

//direct integration of the Lorentz force, yields non-physical energy gain
void updateVelocityDirect(Particle &part, double dt, const double3 &E, const double3 &B) 
{
	double qm = part.q/part.m;		//q/m
	double3 vxB = cross(part.v,B);
	part.v += qm*dt*(E+vxB);
}

//Boris method, conserves energy
void updateVelocityBoris(Particle &part, double dt, const double3 &E, const double3 &B)
{
	double qm = part.q/part.m;		//q/m
    //v minus, first half acceleration
	double3 v_minus = part.v + qm*E*dt/2;

	//v prime, first half rotation
	double3 t = qm*B*dt/2;
	double3 v_prime = v_minus + cross(v_minus, t);

	//v_plus, second half rotation
	double3 s = 2.0*t/(1+dot(t,t));
	double3 v_plus = v_minus + cross(v_prime, s);

	//v[k+0.5], second half acceleration
	part.v = v_plus + qm*E*dt/2;
}

int main() {

	double3 E{0,0,0};		//no electric field
	double3 B{0,0,0.01};	//0.01T in z direction

	//compute Larmor radius for an electron
	double q = -Const::QE;
	double m = Const::ME;
	double vt = 1e3;
	double rL = m*vt/(abs(q)*mag(B));
	
	//create two electrons with orbit centered on origin
	Particle part1({rL,0,0},{0,vt,0},q,m);
	Particle part2({rL,0,0},{0,vt,0},q,m);
	
	double dt = 5e-11;
	
	//output file
	ofstream out("trace.csv");
	out<<"ts p1x p1y p1z p2x p2y p2z\n";

	//rewind velocity per leapfrog
	updateVelocityDirect(part1,-0.5*dt,E,B);
	updateVelocityBoris(part2,-0.5*dt,E,B);
	
	//main loop
	for (int ts = 0;ts<=300;ts++) {		
		//output data
		out<<ts<<" "<<part1.x<<" "<<part2.x<<"\n";		

		updateVelocityDirect(part1,dt,E,B);
		updateVelocityBoris(part2,dt,E,B);
		part1.x += part1.v*dt;
		part2.x += part2.v*dt;	
	}	

return 0;
} 

