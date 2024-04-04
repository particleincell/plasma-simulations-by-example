//Main.cpp
#include <iostream>
#include "World.h"
#include "Output.h"
#include "PotentialSolver.h"
#include "Species.h"

using namespace std;
using namespace Const;

int main() {
  /*initialize domain*/
  	World world(21,21,21);
	world.setExtents({-0.1,-0.1,0.0},{0.1,0.1,0.2});
  
	/*initialize potential solver and solve initial potential*/
    PotentialSolver solver(world,10000,1e-4);
    solver.solve();

    /*obtain initial electric field*/
    solver.computeEF();

	/*set up particle species*/
	vector<Species> species;
	species.push_back(Species("O+", 16*AMU, QE, world));
	species.push_back(Species("e-", ME, -1*QE, world));
	
	/*load particles*/
	
	/*
	int np_ions = 80000;			//number of simulation ions
	int np_eles = 10000;			//number of simulation electrons
	species[0].loadParticlesBox(world.getX0(),world.getXm(),1e11,np_ions);	//ions
	species[1].loadParticlesBox(world.getX0(),world.getXc(),1e11,np_eles);	//electrons
	*/

	int3 np_ions_grid = {21,21,21};
	int3 np_eles_grid = {41,41,41};
	species[0].loadParticlesBoxQS(world.getX0(),world.getXm(),1e11,np_ions_grid);	//ions
	species[1].loadParticlesBoxQS(world.getX0(),world.getXc(),1e11,np_eles_grid);	//electrons

	/*use for-each loop to iterate over all members of a vector*/
	for (Species &sp:species) 
		cout<<sp.name<<" has "<<sp.getNp()<<" particle"<<endl;
	
	/*compute number density*/
	for (Species &sp:species)
		sp.computeNumberDensity();

	/*compute charge density*/
	world.computeChargeDensity(species);

  	Output::fields(world, species);
  	return 0;
}


