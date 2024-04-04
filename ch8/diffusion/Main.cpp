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
#include "Collisions.h"

using namespace std;		//to avoid having to write std::cout
using namespace Const;		//to avoid having to write Const::ME

/*program execution starts here*/
int main(int argc, char *args[])
{
	/*initialize domain*/
    World world(21,21,41);
    world.setExtents({-0.1,-0.1,0},{0.1,0.1,0.4});
    world.setTime(1e-7,10000);

	/*set objects*/
	double phi_sphere = -100;		//set default
	if (argc>1)
		phi_sphere = atof(args[1]);	//convert argument to float
	cout<<"Sphere potential: "<<phi_sphere<<" V"<<endl;
    world.addSphere({0,0,0.15},0.05,phi_sphere);
    world.addInlet();

	/*set up particle species*/
    vector<unique_ptr<Species>> species;
    species.emplace_back(new KineticSpecies("O+", 16*AMU, QE, 5e3, world));
	species.emplace_back(new FluidSpecies("Sph", 100*AMU, 0, 1e16, 100, world));
	Species &ions = *species[0];

	/*setup injection sources*/
	const double ndi = 1e10;			//mean ion density
	vector<unique_ptr<Source>> sources;
	sources.emplace_back(new ColdBeamSource(ions,world,7000,ndi));	//ion source

	/*setup material interactions*/
	vector<unique_ptr<Interaction>> interactions;
	//interactions.emplace_back(new MCC_CEX(ions,neutrals,world));

	/*initialize potential solver and solve initial potential*/
    PotentialSolver solver(world,SolverType::PCG,1000,1e-4);
    solver.setReferenceValues(0,1.5,ndi);
    solver.solve();

    /*obtain initial electric field*/
    solver.computeEF();

    /* main loop*/
	while(world.advanceTime())
    {
		/*inject particles*/
    	for (auto &source:sources)
    		source->sample();

    	/*perform material interactions*/
    	for (auto &interaction:interactions)  interaction->apply(world.getDt());

		/*move particles*/
		for (auto &sp:species)
		{
			sp->advance();
		}

		/*compute charge density*/
		world.computeChargeDensity(species);

        /*update potential*/
       // solver.solve();

        /*obtain electric field*/
        solver.computeEF();

        /*update averages at steady state*/
        if (world.steadyState(species)) {
        	for (auto &sp:species)
        		sp->updateAverages();
        }

		/*screen and file output*/
        Output::screenOutput(world,species);

		/*periodically write out results*/
        if (world.getTs()%50==0 || world.isLastTimeStep()) {
			Output::fields(world, species);
        }

    }
	
	/* grab starting time*/
	cout<<"Simulation took "<<world.getWallTime()<<" seconds"<<endl;
	return 0;		//indicate normal exit
}
