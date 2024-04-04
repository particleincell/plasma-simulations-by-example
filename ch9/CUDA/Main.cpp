#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <thread>
#include "World.h"
#include "PotentialSolver.h"
#include "Species.h"
#include "Output.h"
#include "Source.h"

using namespace std;		//to avoid having to write std::cout
using namespace Const;		//to avoid having to write Const::ME

/*program execution starts here*/
int main(int argc, char *args[])
{
	/*initialize domain*/
    World world(21,21,41);
    world.setExtents({-0.1,-0.1,0},{0.1,0.1,0.4});
    world.setTime(1e-7,400);

	/*check command line arguments for thread count*/
	int num_threads = thread::hardware_concurrency();
	if (argc>1) num_threads = atoi(args[1]);
	cout<<"Running with "<<num_threads<<" threads"<<endl;
	world.setNumThreads(num_threads);   //set number of threads to use

	/*set objects*/
	double phi_sphere = -100;		//set default
	cout<<"Sphere potential: "<<phi_sphere<<" V"<<endl;
    world.addSphere({0,0,0.15},0.05,phi_sphere);
    world.addInlet();

	/*set up particle species*/
	vector<Species> species;
	species.push_back(Species("O", 16*AMU, 0, 1e5, world));
	species.push_back(Species("O+", 16*AMU, QE, 5e0, world));
	species.push_back(Species("O++", 16*AMU, 2*QE, 5e0, world));
	Species &neutrals = species[0];
	Species &ions = species[1];
	Species &ions2 = species[2];
	
	/*setup injection sources*/
	const double nda = 1e13;			//neutral density
	const double ndi = 1e10;			//mean ion density
	vector<ColdBeamSource> sources;
	sources.emplace_back(neutrals,world,7000,nda);		//neutral source
	sources.emplace_back(ions,world,7000,0.8*ndi);	//ion source
	sources.emplace_back(ions2,world,7000,0.1*ndi);	//ion++ source

	/*initialize potential solver and solve initial potential*/
    PotentialSolver solver(world,SolverType::QN,10000,1e-4);
    solver.setReferenceValues(0,1.5,5e12);
    solver.solve();

    /*obtain initial electric field*/
    solver.computeEF();

    /* main loop*/
	while(world.advanceTime())
    {
	  	/*inject particles*/
    	for (ColdBeamSource &source:sources) 
    		source.sample();

        /*move particles*/
		for (Species &sp:species)
		{
			sp.advance();
			sp.computeNumberDensity();
		}

		/*compute charge density*/
		world.computeChargeDensity(species);

        /*update potential*/
        solver.solve();

        /*obtain electric field*/
        solver.computeEF();

        /*update averages at steady state*/
        if (world.steadyState(species)) {
        	for (Species &sp:species)
        		sp.updateAverages();
        }

		/*screen and file output*/
        Output::screenOutput(world,species);
        Output::diagOutput(world,species);

		/*periodically write out results*/
        if (world.getTs()%20==0 || world.isLastTimeStep()) {
        	solver.updateHostPhi();
			Output::fields(world, species);
        }

    }
	
	/* grab starting time*/
	cout<<"Simulation took "<<world.getWallTime()<<" seconds"<<endl;
	return 0;		//indicate normal exit
}
