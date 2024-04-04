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
	double mu = -1;		//
	Particle(const double3 &x, const double3 &v, double q, double m):
	x{x}, v{v}, q{q}, m{m} {}
};

double3 evalB(double z_in) {
    double z = z_in/1e-5;
    double Bz = exp(-2*(z-1)*(z-1)) + exp(-2*(z+1)*(z+1));
    return {0,0,Bz};
}

double3 evalGradB(double z_in) {
    double z = z_in/1e-5;
    double gBz = exp(-2*(z-1)*(z-1))*(-4*(z-1)) + exp(-2*(z+1)*(z+1))*(-4*(z+1));
    return {0,0,gBz*1e5};  //because of unit conversion
}

//Boris method
double updateVelocityBoris(Particle &part, double dt, const double3 &E)
{
	double qm = part.q/part.m;		//q/m

	//evaluate B field and gradB at particle position
	double3 B = evalB(part.x[2]);
	double3 gradB = evalGradB(part.x[2]);

	//unit vector in B direction
	double3 B_unit = unit(B);
	

	double mag1, mag2;
	//magnetic moment, should be conserved
	if (part.mu<0) {
		double3 v_par = dot(part.v,B_unit)*B_unit;
		double3 v_perp = part.v - v_par;	
		part.mu=0.5*part.m*dot(v_perp,v_perp)/mag(B);	//initialize
	}		
	else {
	//orbital velocity, part.v - v_parallel	
		mag1 = mag(part.v);	
		
		double3 v_par = dot(part.v,B_unit)*B_unit;
		double3 v_perp = part.v - v_par;
		double v_perp_mag_new = sqrt(part.mu/part.m*2*mag(B));
		double v_perp_mag = mag(v_perp);
		v_perp = v_perp_mag_new*unit(v_perp);
		double v_par2 =mag1*mag1 - v_perp_mag_new*v_perp_mag_new;
		double v_par_mag_new;		
		if (v_par2>=0) v_par_mag_new = sqrt(v_par2);
		else v_par_mag_new = -sqrt(-v_par2);
		v_par = v_par_mag_new * unit(v_par);
		part.v = v_perp+v_par;	
		mag2 = mag(part.v);

	}

	//total acceleration force
	double3 F = part.q*E; //- mu*gradB*unit(gradB);
	
    //v minus, first half acceleration
	double3 v_minus = part.v + (F/part.m)*dt/2;

	//v prime, first half rotation
	double3 t = qm*B*dt/2;
	double3 v_prime = v_minus + cross(v_minus, t);

	//v_plus, second half rotation
	double3 s = 2.0*t/(1+dot(t,t));
	double3 v_plus = v_minus + cross(v_prime, s);

	//v[k+0.5], second half acceleration
	part.v = v_plus + (F/part.m)*dt/2;
//cout<<part.x[2]<<": "<<sqrt(part.v[0]*part.v[0]+part.v[1]*part.v[1]);

//cout<<" ZZ "<<part.v[2]<<" ("<<mag1<<", "<<mag2<<")"<<endl;
		double3 v_par = dot(part.v,B_unit)*B_unit;
		double3 v_perp = part.v - v_par;
	return 0.5*part.m*dot(v_perp,v_perp)/mag(B);
}

int main() {

	double3 E{0,0,0};		//no electric field

	//compute Larmor radius for an electron
	double q = -Const::QE;
	double m = Const::ME;
	double vt = 1e4;
	double3 B = evalB(0);
	double rL = m*vt/(abs(q)*mag(B));
	
	//create two electrons with orbit centered on origin
	Particle part({rL,0,0},{0,vt,1.6e4},q,m);
	double dt = 1e-12;
	
	//output file
	ofstream out("trace.csv");
	out<<"ts x y z u v w v_perp v_par mag(v) Bz gradBz mu\n";

	double mu;

	//rewind velocity per leapfrog
	mu = updateVelocityBoris(part,-0.5*dt,E);
	
	//main loop
	for (int ts = 0;ts<=20000;ts++) {		
		//output data
		out<<ts<<" "<<part.x<<" "<<part.v<<" "<<sqrt(part.v[0]*part.v[0]+part.v[1]*part.v[1])<<" "<<part.v[2]<<" " <<mag(part.v)<<" "<<evalB(part.x[2])[2]<<" "<<evalGradB(part.x[2])[2]<<" "<<mu<<"\n";		

		mu = updateVelocityBoris(part,dt,E);
		part.x += part.v*dt;
	}	

return 0;
} 

