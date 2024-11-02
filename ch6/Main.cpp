#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include "World.h"
#include "FESolver.h"
#include "Species.h"
#include "Output.h"
#include "Source.h"

using namespace std;		//to avoid having to write std::cout
using namespace Const;		//to avoid having to write Const::ME

/*program execution starts here*/
int main(int argc, char *args[])
{
	/*initialize domain*/
	/*instantiate volume*/
	VolumeMesh vm("mesh.unv");

	World world(vm);
	world.setTime(1e-7,800);

	/*set up particle species*/
	vector<Species> species;
	species.emplace_back("O+", 16*AMU, QE, 1e2, world);

	/*setup injection sources*/
	vector<ColdBeamSource> sources;
	Group &group = world.vm.groups.at("inlet"); //throws an exception if not found
	sources.emplace_back(species[0],world,7000,1e10,group);	//ion source

	

	/*initialize potential solver and solve initial potential*/
    FESolver solver(world,10000,1e-4);
    solver.setReferenceValues(0,1.5,1e10);
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
			sp.advance(solver);
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
        {
			Output::fields(world, species);
			Output::particles(world, species, 10000);
        }
    }

	/* grab starting time*/
	cout<<"Simulation took "<<world.getWallTime()<<" seconds"<<endl;
	return 0;		//indicate normal exit
}
