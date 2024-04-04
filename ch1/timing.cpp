/*
Force calculation timing
*/

#include <iostream>
#include <cstring>
#include <sstream>
#include <chrono>
#include <ctime>
#include <math.h>
#include <random>

/*constants*/
#define EPS0 8.85418782e-12  	// F/m, vacuum permittivity
#define K	1.38065e-23			// J/K, Boltzmann constant
#define ME 9.10938215e-31		// kg, electron mass
#define QE 1.602176565e-19		// C, electron charge
#define AMU  1.660538921e-27	// kg, atomic mass unit


//to avoid having to write std::cout and so on
using namespace std;
using namespace std::chrono;

std::mt19937 mt_gen(0);		/*seed*/
std::uniform_real_distribution<double> rnd_dist(0, 1.0);
double rnd()
{
	return rnd_dist(mt_gen);
}

/* Data structure for particle storage **/
struct Particle
{
	double x[3];			/*position*/	
	Particle()
	{
		x[0] = rnd();
		x[1] = rnd();
		x[2] = rnd();
	}
};


/* --------- main -------------*/
int main()
{
	int np = 100000000;
	Particle *particles = new Particle[np];
	
	//grab starting time clock
	using namespace std::chrono;
	high_resolution_clock::time_point t_start = high_resolution_clock::now();
	/*note we could also simplify this by using "auto" variable as in
	  auto t_start = high_resolution_clock::now();	*/
	  
	 double F[3] = {0,0,0};
	 double r_vec[3];
	 
	 Particle *self = &particles[0];
	 double pi = acos(-1);
	 double kc = 1/(4*pi*EPS0);
	 
	 for (Particle *part = &particles[1]; part<&particles[np];part++)
	 {
	 	r_vec[0] = part->x[0]-self->x[0];
	 	r_vec[1] = part->x[1]-self->x[1];
	 	r_vec[2] = part->x[2]-self->x[2];
		double r_mag = sqrt(r_vec[0]*r_vec[0]+r_vec[1]*r_vec[1]+r_vec[2]*r_vec[2]);
		double r_mag3 = r_mag*r_mag*r_mag;
		F[0] += kc*(QE*QE)*r_vec[0]/r_mag3;
		F[1] += kc*(QE*QE)*r_vec[0]/r_mag3;
		F[2] += kc*(QE*QE)*r_vec[0]/r_mag3;		
	 }
	 
	//grab ending time
	high_resolution_clock::time_point t_end = high_resolution_clock::now();
    std::chrono::duration<double> duration = t_end-t_start;
	cout<<"Simulation took "<<duration.count()/(np-1)<<" s per particle"<<endl;
	
	//so the loop doesn't get optimized away
	cout<<"Total Force: "<<F[0]<<", "<<F[1]<<", "<<F[2]<<endl;
	
	return 0;
}
