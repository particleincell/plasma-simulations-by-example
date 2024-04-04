/*defines the simulation domain*/
#include <random>
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include "World.h"
#include "Field.h"
#include "Species.h"
#include "FESolver.h"
	
//make an instance of the Rnd class
Rnd rnd;

using namespace std;

/*constructor*/
World::World(VolumeMesh &vm): vm{vm}, node_vol(vm), phi(vm),rho(vm), ef(vm) {
	time_start =  chrono::high_resolution_clock::now();	//save starting time point

	/*compute node volumes by summing up fractional cell volumes*/
	for (Tet &tet: vm.tets) {
		for (int v=0;v<4;v++)
			node_vol[tet.con[v]] += 0.25*tet.volume;
	}
}

/*returns elapsed wall time in seconds*/
double World::getWallTime() {
  auto time_now = chrono::high_resolution_clock::now();
  chrono::duration<double> time_delta = time_now-time_start;
  return time_delta.count();
}

/*computes charge density from rho = sum(charge*den)*/
void World::computeChargeDensity(vector<Species> &species)
{
	rho = 0;
	for (Species &sp:species)
	{
		if (sp.charge==0) continue;	//don't bother with neutrals
		rho += sp.charge*sp.den;
	}
}

/* computes total potential energy from 0.5*eps0*sum(E^2)*/
double World::getPE() {
	double pe = 0;
	for (size_t i=0;i<vm.nodes.size();i++)
	{
		double3 efn = ef[i];	//ef at this node
		double ef2 = efn[0]*efn[0]+efn[1]*efn[1]+efn[2]*efn[2];

		pe += ef2*node_vol[i];
	}
	return 0.5*Const::EPS_0*pe;
}

/*brute force search*/
LCord World::XtoLbrute(const double3 &pos) {
	LCord lc;
	for (size_t e = 0;e<vm.tets.size();e++)
	{
		LCord lc = vm.XtoL(pos,e,false);
		if (lc) return lc;
	}
	return LCord(-1);
}

/*checks for steady state by comparing change in mass, momentum, and energy*/
bool World::steadyState(vector<Species> &species) {
	// do not do anything if already at steady state
	if (steady_state) return true;

	double tot_mass = 0;
	double tot_mom = 0;
	double tot_en = getPE();
	for (Species &sp:species)
	{
		tot_mass += sp.getRealCount();	//number of real molecules
		double3 mom = sp.getMomentum();
		tot_mom += mag(mom);		//z-component of momentum
		tot_en += sp.getKE();		//add kinetic energy
	}

	/*compute new values to last*/
	const double tol = 1e-3;
	if (abs((tot_mass-last_mass)/tot_mass)<tol &&
		abs((tot_mom-last_mom)/tot_mom)<tol &&
		abs((tot_en-last_en)/tot_en)<tol) {
		steady_state = true;
		cout<<"Steady state reached at time step "<<ts<<endl;
	}

	/*update prior values*/
	last_mass = tot_mass;
	last_mom = tot_mom;
	last_en = tot_en;
	return steady_state;
}



