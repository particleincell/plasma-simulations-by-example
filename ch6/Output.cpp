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
	VolumeMesh &vm = world.vm;

	stringstream ss;
	ss<<"results/mesh_"<<setfill('0')<<setw(4)<<world.getTs()<<".vtu";
	ofstream out(ss.str());
	if (!out.is_open()) {cerr<<"Failed to open file "<<ss.str()<<endl;exit(-1);}

	/*header*/
	out<<"<?xml version=\"1.0\"?>\n";
	out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out<<"<UnstructuredGrid>\n";
	out<<"<Piece NumberOfPoints=\""<<vm.nodes.size()<<"\" NumberOfVerts=\"0\" NumberOfLines=\"0\" ";
	out<<"NumberOfStrips=\"0\" NumberOfCells=\""<<vm.tets.size()<<"\">\n";
	
	/*points*/
	out<<"<Points>\n";
	out<<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for (Node &node: vm.nodes)
		out<<node.pos[0]<<" "<<node.pos[1]<<" "<<node.pos[2]<<"\n";
	out<<"</DataArray>\n";
	out<<"</Points>\n";

	/*Cells*/
	out<<"<Cells>\n";
	out<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
	for (Tet &tet: vm.tets)
		out<<tet.con[0]<<" "<<tet.con[1]<<" "<<tet.con[2]<<" "<<tet.con[3]<<"\n";
	out<<"</DataArray>\n";

	out<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	for (size_t e=0; e<vm.tets.size();e++)
		out<<(e+1)*4<<" ";
	out<<"\n";
	out<<"</DataArray>\n";

	out<<"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
	for (size_t e=0; e<vm.tets.size();e++)
		out<<"10 ";		//tetrahedron
	out<<"\n";
	out<<"</DataArray>\n";
	out<<"</Cells>\n";

	/*save point data*/
	out<<"<PointData Scalars=\"phi\">\n";

	//node index
	out<<"<DataArray type=\"Int32\" Name=\"node_index\" format=\"ascii\">\n";
	for (size_t n=0; n<vm.tets.size();n++)	out<<n<<" ";
	out<<"\n</DataArray>\n";

	//node types
	out<<"<DataArray type=\"Int32\" Name=\"node_type\" format=\"ascii\">\n";
	for (size_t n=0; n<vm.nodes.size();n++)	out<<vm.nodes[n].type<<" ";
	out<<"\n</DataArray>\n";

	//node volumes
	out<<"<DataArray type=\"Float64\" Name=\"node_vol\" format=\"ascii\">\n";
	out<<world.node_vol;
	out<<"</DataArray>\n";


	out<<"<DataArray type=\"Float64\" Name=\"phi\" format=\"ascii\">\n";
	out<<world.phi;
	out<<"</DataArray>\n";
	
	/*output species densities*/
	for (Species &sp:species)
	{
		out<<"<DataArray Name=\"nd."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.den;
		out<<"</DataArray>\n";

		out<<"<DataArray Name=\"nd-ave."<<sp.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		out<<sp.den_ave;
		out<<"</DataArray>\n";

	}

	out<<"</PointData>\n";

	/*save cell data*/
	out<<"<CellData Vectors=\"ef\">\n";

	out<<"<DataArray type=\"Float64\" Name=\"cell_volume\" format=\"ascii\">\n";
	for (Tet &tet:world.vm.tets) out<<tet.volume<<" ";
	out<<"\n</DataArray>\n";

	out<<"</CellData>\n";

	out<<"</Piece>\n";
	out<<"</UnstructuredGrid>\n";
	out<<"</VTKFile>\n";

	out.close();
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


/*saves particle data*/
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
				{to_output.emplace_back(&part);counter=0;}
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
