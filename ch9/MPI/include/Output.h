/*output in VTK format*/
#ifndef _OUTPUT_H
#define _OUTPUT_H

#include "World.h"
#include "Species.h"
#include <string>
#include <vector>

using namespace std;

class Output 
{
public:
 	
	Output(World &world) : world(world) {}

    void close()
	{
		log.close();
    }
    
    /*saves simulation state data vs time*/
    void outputLog(int ts, int num_colls, vector<Species> &species_list);
	
    /*saves field data*/
    void outputField(int ts, vector<Species> &species_list);
    
    /*outputs trace data - the same particle is saved at each time*/
    void outputTracer(int ts, Species &species, int count);
    
    /*outputs scatter data, random sample of particles*/
    void outputScatter(int ts, Species &species, int count);

    /*saves pvd file to animate multiblock data*/
    void outputPVD();

/*data*/
	World &world;
    string path="results/";	
    ofstream log;
    int last_field_output=-99;

    vector<int>save_times;	/*array of time steps when data was saved, used to make pvd*/
};

#endif
