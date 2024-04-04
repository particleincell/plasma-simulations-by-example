#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <memory>
#include "World.h"
#include "EMSolver.h"
#include "MagneticSolver.h"
#include "Species.h"
#include "Output.h"
#include "Source.h"

using namespace std;		//to avoid having to write std::cout
using namespace Const;		//to avoid having to write Const::ME

/*program execution starts here*/
int main(int argc, char *args[])
{
   /*initialize domain*/
    World world(41,41);
    world.setExtents({0,-0.02},{0.04,0.02});
    double dt = world.getDh()[0]/Const::C;		//dx = 0.5*c*dt
    world.setTime(dt,10000);
    cout<<"dt = "<<dt<<endl;

	/*set objects*/
	double phi_circle = -0.001;		//set default
	world.addCircle({0.02,0},0.0015,phi_circle);
    world.addInlet();

	/*set up particle species*/
    vector<Species> species;
    species.push_back(Species("e-", ME, -1*QE, 1e-0, world));
    species.push_back(Species("H+", AMU, 1*QE, 1e-0, world));
    Species &eles = species[0];
	Species &ions = species[1];
	
	/*setup injection sources*/
	const double nde = 2e9;			//mean ion density
	vector<ColdBeamSource> sources;
	sources.emplace_back(eles,world,5e6,nde);	//electron source
	sources.emplace_back(ions,world,1e6,nde);	//ion source

    /*magnetic field solution*/
    //world.magnetizeSphere({0,100000,0});
    //MagneticSolver bsolver(world,20000,1e-4);
    //bsolver.solve();

	/*initialize potential solver and solve initial potential*/
    EMSolver solver(world,10000,1e-3);
    solver.advance();

	/* main loop*/
	while(world.advanceTime())
    {
		/*inject particles*/
    	for (auto source:sources) {
    		source.sample();
    	}

		/*move particles*/
    	for (Species &sp:species)
		{
			sp.advance();
			sp.computeGasProperties();
		}

		/*compute charge density*/
		world.computeChargeDensity(species);
		world.computeCurrentDensity(species);

        /*update fields*/
        solver.advance();

        /*update averages at steady state*/
        if (world.steadyState(species)) {
        	for (Species &sp:species)
        		sp.updateAverages();
        }

		/*screen and file output*/
        Output::screenOutput(world,species);
        Output::diagOutput(world,species);

		/*periodically write out results*/
        if (world.getTs()%100==0 || world.isLastTimeStep())
			Output::fields(world, species);
    }
	
	/* show run time*/
	cout<<"Simulation took "<<world.getWallTime()<<" seconds"<<endl;
	return 0;		//indicate normal exit
}
