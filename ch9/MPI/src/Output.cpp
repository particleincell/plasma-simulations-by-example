#include "Output.h"
#include <list>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;

/*saves simulation state (number of particles, etc..) in csv format*/
void Output::outputLog(int ts, int num_colls, vector<Species> &species_list)
{
	/*open log file if not yet open*/
    if (!log.is_open()) 
	{
		log.open("log.cvs");
		if (!log.is_open()) {cerr<<"Failed to open log file!"<<endl; return;}

		/*write header*/
		log<<"ts,time";
        for (Species &species:species_list)
                log<<",np."<<species.name<<",mom."<<species.name;
        log<<",num_cols\n";
	}

    log<<ts<<","<<world.time;
    for (Species &species:species_list)
        log<<","<<species.np<<","<<species.getMomentum();
    log<<","<<num_colls<<"\n";    

	if (ts%100==0) log.flush();
}

/*saves output in VTK format*/
void Output::outputField(int ts, vector<Species> &species_list)
{
	/*do not output if already save at this time step*/
    if (ts==last_field_output) return;
    last_field_output = ts;

	stringstream name;
	name<<"results_"<<setfill('0')<<setw(5)<<ts<<"_"<<world.mpi_rank<<".vti";

    /*open output file*/
    ofstream out(path+name.str());
   	if (!out.is_open()) {cerr<<"Could not open "<<name.str()<<endl;return;}

	/*ImageData is vtk format for structured Cartesian meshes*/
	out<<"<VTKFile type=\"ImageData\">\n";
	out<<"<ImageData Origin=\""<<world.x0[0]<<" "<<world.x0[1]<<" "<<world.x0[2]<<"\" ";
	out<<"Spacing=\""<<world.dh[0]<<" "<<world.dh[1]<<" "<<world.dh[2]<<"\" ";
	out<<"WholeExtent=\"0 "<<world.ni-1<<" 0 "<<world.nj-1<<" 0 "<<world.nk-1<<"\">\n";
	
	/*output data stored on nodes (point data)*/
	out<<"<PointData>\n";
	
	/*potential, scalar*/
	out<<"<DataArray Name=\"phi\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	world.phi->write(out);
	out<<"</DataArray>\n";

	/*electric field, 3 component vector*/
	out<<"<DataArray Name=\"ef\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
	world.ef->write(out);
	out<<"</DataArray>\n";
	
	/*output species densities*/
	for (Species &species:species_list)
	{
		out<<"<DataArray Name=\"nd."<<species.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		species.den->write(out);
		out<<"</DataArray>\n";
	}

	/*repeat with average values*/
	for (Species &species:species_list)
	{
		out<<"<DataArray Name=\"nd-ave."<<species.name<<"\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
		if (species.den_ave->num_samples>0)		/*if we have averaged data*/	
			species.den_ave->write(out);
		else
			species.den->write(out);
		out<<"</DataArray>\n";
	}

	/*output species velocities*/
	for (Species &species:species_list)
	{
		out<<"<DataArray Name=\"vel-ave."<<species.name<<"\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
		if (species.vel_ave.f[0]->num_samples>0)	/*if we have averaged data*/
			species.vel_ave.write(out);
		else
			species.vel.write(out);
		out<<"</DataArray>\n";
	}
		
	/*close out tags*/
	out<<"</PointData>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";
 	out.close();

 	/*record this time step*/
 	save_times.emplace_back(ts);
}

void Output::outputPVD()
{
	/*also output pvd file on rank 0*/
 	if (world.mpi_rank!=0) return;

 	stringstream name;
	ofstream out("results.pvd");
	out<<"<?xml version=\"1.0\"?>\n";
	out<<"<VTKFile type=\"Collection\" version=\"0.1\">\n";
	out<<"<Collection>\n";
	for (int t=0;t<(int)save_times.size();t++)
	{
		int ts = save_times[t];
		for (int i=0;i<world.mpi_size;i++)
		{
			stringstream name;
			name<<path<<"results_"<<setfill('0')<<setw(5)<<ts<<"_"<<i<<".vti";
			out<<"<DataSet timestep=\""<<t<<"\" file=\""<<name.str()<<"\"/>\n";
		}
	}
	out<<"</Collection>\n";
	out<<"</VTKFile>\n";
	out.close();
}

/*outputs scatter data, random sample of particles*/
void Output::outputScatter(int ts, Species &species, int count)
{
    stringstream name;
	name<<"scatter_"<<setfill('0')<<setw(5)<<ts<<".vts";

    /*open output file*/
    ofstream out(path+name.str());
   	if (!out.is_open()) {cerr<<"Could not open "<<name.str()<<endl;return;}
        
        
    /*generate the random sample*/
	double fraction = count/(double)species.np;	//fraction of particles to save
    
	vector<Particle*> picks;
				
	for (int p=0;p<species.np;p++)
    {
		if (fraction>=rnd())
			picks.push_back(&species.particles[p]);
    }
	
    /*StructuredGrid has logical topology but arbitrarily located nodes*/
	out<<"<VTKFile type=\"StructuredGrid\">\n";
	out<<"<StructuredGrid WholeExtent=\"0 "<<picks.size()-1<<" 0 0 0 0\">\n";
	
	/*output data stored on nodes (point data)*/
	out<<"<PointData>\n";

	/*particle ids*/
	out<<"<DataArray Name=\"id\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	for (Particle *part:picks) out<<part->id<<" ";
	out<<"\n</DataArray>\n";

	/*velocities*/
	out<<"<DataArray Name=\"vel\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
	for (Particle *part:picks) out<<part->vel[0]<<" "<<part->vel[1]<<" "<<part->vel[2]<<"\n";	
	out<<"</DataArray>\n";
	out<<"</PointData>\n";

	/*"node" positions*/
	out<<"<Points>\n";
	out<<"<DataArray Name=\"pos\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
	for (Particle *part:picks) out<<part->pos[0]<<" "<<part->pos[1]<<" "<<part->pos[2]<<"\n";
	out<<"</DataArray>\n";
	out<<"</Points>\n";

	out<<"</StructuredGrid>\n";
	out<<"</VTKFile>\n";
	out.close();
}
    
/*like scatter, but outputs first "count" particles*/
void Output::outputTracer(int ts, Species &species, int count)
{
    stringstream name;
	name<<"tracer_"<<setfill('0')<<setw(5)<<ts<<".vts";

    /*open output file*/
    ofstream out(path+name.str());
   	if (!out.is_open()) {cerr<<"Could not open "<<name.str()<<endl;return;}
        
    vector<Particle*> picks;
				
    /*make sure count is not more than particle count*/
	int min_count = species.np;
	if (count<min_count) min_count = count;

	for (int i=0;i<min_count;i++)
    {
		picks.push_back(&species.particles[i]);
    }
	
    /*StructuredGrid has logical topology but arbitrarily located nodes*/
	out<<"<VTKFile type=\"StructuredGrid\">\n";
	out<<"<StructuredGrid WholeExtent=\"0 "<<picks.size()-1<<" 0 0 0 0\">\n";
	
	/*output data stored on nodes (point data)*/
	out<<"<PointData>\n";

	/*particle ids*/
	out<<"<DataArray Name=\"id\" NumberOfComponents=\"1\" format=\"ascii\" type=\"Float64\">\n";
	for (Particle *part:picks) out<<part->id<<" ";
	out<<"\n</DataArray>\n";

	/*velocities*/
	out<<"<DataArray Name=\"vel\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
	for (Particle *part:picks) out<<part->vel[0]<<" "<<part->vel[1]<<" "<<part->vel[2]<<"\n";	
	out<<"</DataArray>\n";
	out<<"</PointData>\n";

	/*"node" positions*/
	out<<"<Points>\n";
	out<<"<DataArray Name=\"pos\" NumberOfComponents=\"3\" format=\"ascii\" type=\"Float64\">\n";
	for (Particle *part:picks) out<<part->pos[0]<<" "<<part->pos[1]<<" "<<part->pos[2]<<"\n";
	out<<"</DataArray>\n";
	out<<"</Points>\n";
	
	out<<"</StructuredGrid>\n";
	out<<"</VTKFile>\n";
	out.close();
}
    

