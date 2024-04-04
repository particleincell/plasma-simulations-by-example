#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
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

	/*set objects*/
	double phi_sphere = -100;		//set default
	if (argc>1)
		phi_sphere = atof(args[1]);	//convert argument to float
	cout<<"Sphere potential: "<<phi_sphere<<" V"<<endl;
    world.addSphere({0,0,0.15},0.05,phi_sphere);
    world.addInlet();

	/*set up particle species*/
	vector<Species> species;
	species.push_back(Species("O+", 16*AMU, QE, 1e4, world));

	/*setup injection sources*/
	vector<ColdBeamSource> sources;
	sources.push_back(ColdBeamSource(species[0],world,7000,1e12));	//ion source

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
        if (world.getTs()%20==0 || world.isLastTimeStep())
			Output::fields(world, species);

    }
	
	/* grab starting time*/
	cout<<"Simulation took "<<world.getWallTime()<<" seconds"<<endl;
	return 0;		//indicate normal exit
}
