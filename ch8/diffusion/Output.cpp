#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "Output.h"
#include "World.h"
#include "Species.h"

using namespace std;

/*saves fields in VTK format*/
void Output::fields(World &world, vector<unique_ptr<Species>>  &species)
{
	/*update gas macroscopic properties*/
	for (auto &sp:species)
		sp->computeGasProperties();

	stringstream name;
	name<<"results/fields_"<<setfill('0')<<setw(5)<<world.getTs()<<".vti";

    /*open output file*/
    ofstream out(name.str());
   	if (!out.is_open()) {cerr<<"Could not open "<<name.str()<<endl;return;}

	/*ImageData is vtk format for structured Cartesian meshes*/
	out<<"<VTKFile type=\"ImageData\">\n";
	double3 x0 = world.getX0();
	double3 dh = world.getDh();
	out<<"<ImageData Origin=\""<<x0[0]<<" "<<x0[1]<<" "<<x0[2]<<"\" ";
	out<<"Spacing=\""<<dh[0]<<" "<<dh[1]<<" "<<dh[2]<<"\" ";
	out<<"WholeExtent=\"0 "<<world.ni-1<<" 0 "<<world.nj-1<<" 0 "<<world.nk-1<<"\">\n";
	
	/*output data stored on nodes (point data)*/
	out<<"<PointData>\n";

	/*object id, scalar*/
	out<<"<DataArray Name=\"object_id\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Int32\">\n";
	out<<world.object_id;
	out<<"</DataArray>\n";

	/*node volumes, scalar*/
	Field f(world.node_vol);
	for (int i=0;i<world.ni;i++)
		for (int j=0;j<world.nj;j++)
			for (int k=0;k<world.nk;k++)
				f[i][j][k] = world.node_vol[i][j][k];
	out<<"<DataArray Name=\"NodeVol\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.node_vol;
	out<<"</DataArray>\n";

	/*potential, scalar*/
	out<<"<DataArray Name=\"phi\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.phi;
	out<<"</DataArray>\n";

	/*charge density, scalar*/
	out<<"<DataArray Name=\"rho\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.rho;
	out<<"</DataArray>\n";

	/*species number densities*/
	for (auto &sp:species)
	{
		out<<"<DataArray Name=\"nd."<<sp->name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp->den;
		out<<"</DataArray>\n";
	}
	
	/*species averaged number densities*/
	for (auto &sp:species)
	{
		out<<"<DataArray Name=\"nd-ave."<<sp->name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp->den_ave;
		out<<"</DataArray>\n";
	}

	/*species stream velocity*/
	for (auto &sp:species)
	{
		out<<"<DataArray Name=\"vel."<<sp->name<<"\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp->vel;
		out<<"</DataArray>\n";
	}

	/*species temperature*/
	for (auto &sp:species)
	{
		out<<"<DataArray Name=\"T."<<sp->name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp->T;
		out<<"</DataArray>\n";
	}

	/*electric field, 3 component vector*/
	out<<"<DataArray Name=\"ef\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
	out<<world.ef;
	out<<"</DataArray>\n";

	/*close out tags*/
	out<<"</PointData>\n";

	/*cell data*/
	out<<"<CellData>\n";
	/*species temperature*/
	for (auto &sp:species)
	{
		KineticSpecies *ks = dynamic_cast<KineticSpecies*>(sp.get());
		if (!ks) continue; 	//skip non-kinetic species

		out<<"<DataArray Name=\"mpc."<<ks->name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<ks->mpc;
		out<<"</DataArray>\n";
	}
	out<<"</CellData>\n";

	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";
 	out.close();

 	/*clear samples if not at steady state*/
 	if (!world.isSteadyState())
 		for (auto &sp:species) sp->clearSamples();
}

//writes information to the screen
void Output::screenOutput(World &world, vector<unique_ptr<Species>> &species)
{
	cout<<"ts: "<<world.getTs();
	for (auto &sp:species)
		cout<<" "<<sp->printSelf();
	cout<<endl;
}

//file stream handle
namespace Output {
std::ofstream f_diag;
}


