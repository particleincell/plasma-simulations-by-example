#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <memory>
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
    World world(11,11,41);
    world.setExtents({0,0,0},{0.1,0.1,0.4});
    world.setTime(1e-7,400);

	/*set objects*/
	double phi_sphere = -100;		//set default
	world.addSphere({0,0,0.15},0.05,phi_sphere);
    world.addInlet();

	/*set up particle species*/
    vector<Species> species;
    species.push_back(Species("O+", 16*AMU, QE, 5e1, world));
	
	/*setup injection sources*/
	const double ndi = 1e10;			//mean ion density
	vector<ColdBeamSource> sources;
	sources.emplace_back(species[0],world,7000,ndi);	//ion source
	
	/*initialize potential solver and solve initial potential*/
    PotentialSolver solver(world,SolverType::GS,20000,1e-6);
    solver.setReferenceValues(0,1.5,ndi);
    solver.solve();

    /*obtain initial electric field*/
    solver.computeEF();

    /* main loop*/
	while(world.advanceTime())
    {
		/*inject particles*/
    	for (auto source:sources)
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
        if (world.getTs()%50==0 || world.isLastTimeStep()) 
			Output::fields(world, species);
    }
	
	/* show run time*/
	cout<<"Simulation took "<<world.getWallTime()<<" seconds"<<endl;
	return 0;		//indicate normal exit
}
