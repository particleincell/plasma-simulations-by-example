#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "Output.h"
#include "World.h"
#include "Species.h"

using namespace std;

/*saves fields in VTK format*/
void Output::fields(World &world, vector<Species> &species)
{
	/*update gas macroscopic properties*/
	for (Species &sp:species)
		sp.computeGasProperties();

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
	for (Species &sp:species)
	{
		out<<"<DataArray Name=\"nd."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.den;
		out<<"</DataArray>\n";
	}
	
	/*species averaged number densities*/
	for (Species &sp:species)
	{
		out<<"<DataArray Name=\"nd-ave."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.den_ave;
		out<<"</DataArray>\n";
	}

	/*species stream velocity*/
	for (Species &sp:species)
	{
		out<<"<DataArray Name=\"vel."<<sp.name<<"\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.vel;
		out<<"</DataArray>\n";
	}

	/*species temperature*/
	for (Species &sp:species)
	{
		out<<"<DataArray Name=\"T."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.T;
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
	for (Species &sp:species)
	{
		out<<"<DataArray Name=\"mpc."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.mpc;
		out<<"</DataArray>\n";
	}
	out<<"</CellData>\n";

	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";
 	out.close();

 	/*clear samples if not at steady state*/
 	if (!world.isSteadyState())
 		for (Species &sp:species) sp.clearSamples();
}

//writes information to the screen
void Output::screenOutput(World &world, vector<Species> &species)
{
	cout<<"ts: "<<world.getTs();
	for (Species &sp:species)
		cout<<setprecision(3)<<"\t "<<sp.name<<":"<<sp.getNp();
	cout<<endl;
}

//file stream handle
namespace Output {
std::ofstream f_diag;
}

/*save runtime diagnostics to a file*/
void Output::diagOutput(World &world, vector<Species> &species)
{
	using namespace Output;	//to get access to f_diag

	//is the file open?
	if (!f_diag.is_open())
	{
		f_diag.open("runtime_diags.csv");
		f_diag<<"ts,time,wall_time";
		for (Species &sp:species)
			f_diag<<",mp_count."<<sp.name<<",real_count."<<sp.name
				  <<",px."<<sp.name<<",py."<<sp.name<<",pz."<<sp.name
			      <<",KE."<<sp.name;
		f_diag<<",PE,E_total"<<endl;
	}

	f_diag<<world.getTs()<<","<<world.getTime();
	f_diag<<","<<world.getWallTime();

	double tot_KE = 0;
	for (Species &sp:species)
	{
		double KE = sp.getKE();	//species kinetic energy
		tot_KE += KE;		//increment total energy
		double3 mom = sp.getMomentum();

		f_diag<<","<<sp.getNp()<<","<<sp.getRealCount()
			  <<","<<mom[0]<<","<<mom[1]<<","<<mom[2]<<","<<KE;
	}

	//write out system potential and total energy
	double PE = world.getPE();
	f_diag<<","<<PE<<","<<(tot_KE+PE);

	f_diag<<"\n";	//use \n to avoid flush to disc
	if (world.getTs()%25==0) f_diag.flush();
}

/*saves particle x-vel data*/
void Output::particles(World &world, vector<Species> &species, int num_parts) {
	/*loop over all species*/

	for (Species &sp:species) {

		//open a phase_sp_it.vtp
		stringstream name;
		name<<"results/parts_"<<sp.name<<"_"<<setfill('0')<<setw(5)<<world.getTs()<<".vtp";

		/*open output file*/
		ofstream out(name.str());
		if (!out.is_open()) {cerr<<"Could not open "<<name.str()<<endl;return;}

		/*build a list of particles to output*/
		double dp = num_parts/(double)sp.getNp();
		double counter = 0;
		vector<Particle*> to_output;
		for (Particle &part : sp.particles)	{
			counter+=dp;
			if (counter>1) //save particle
				{to_output.emplace_back(&part);counter=-1;}
		}

		/*header*/
		out<<"<?xml version=\"1.0\"?>\n";
		out<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
		out<<"<PolyData>\n";
		out<<"<Piece NumberOfPoints=\""<<to_output.size()<<"\" NumberOfVerts=\"0\" NumberOfLines=\"0\" ";
		out<<"NumberOfStrips=\"0\" NumberOfCells=\"0\">\n";

		/*points*/
		out<<"<Points>\n";
		out<<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
		for (Particle *part: to_output)
			out<<part->pos<<"\n";
		out<<"</DataArray>\n";
		out<<"</Points>\n";

		/*velocities*/
		out<<"<PointData>\n";
		out<<"<DataArray Name=\"vel."<<sp.name<<"\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
		for (Particle *part: to_output)
			out<<part->vel<<"\n";
		out<<"</DataArray>\n";
		out<<"</PointData>\n";

		out<<"</Piece>\n";
		out<<"</PolyData>\n";
		out<<"</VTKFile>\n";

		out.close();
	}

}
