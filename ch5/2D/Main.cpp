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
    World world(41,21);
    world.setExtents({0,-0.1},{0.4,0.1});
    world.setTime(1e-7,400);

	/*set objects*/
	double phi_circle = -100;		//set default
	world.addCircle({0.15,0},0.05,phi_circle);
    world.addInlet();

	/*set up particle species*/
    vector<Species> species;
    species.push_back(Species("O+", 16*AMU, 1*QE, 5e1, world));
	Species &ions = species[0];
	
	/*setup injection sources*/
	const double ndi = 1e10;			//mean ion density
	vector<ColdBeamSource> sources;
	sources.emplace_back(ions,world,7000,ndi);	//ion source
	
	/*initialize potential solver and solve initial potential*/
    PotentialSolver solver(world,SolverType::GS,10000,1e-3);
    solver.setReferenceValues(0,1.5,ndi);
    solver.solve();

    /*obtain initial electric field*/
    solver.computeEF();

    double sigma_cr_max = 1e-16;
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
