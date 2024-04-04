/*Demo 3D Particle In Cell Code

Simulates flow of charged ions past a sphere at a fixed potential

Developed for the Fundamentals of the Particle In Cell Method course
*/

#include "World.h"
#include "Field.h"
#include "Species.h"
#include "PotentialSolver.h"
#include "Source.h"
#include "Output.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <chrono>
#include <thread>
#include <mpi.h>

using namespace std;

/*function prototypes*/
void computeRhoi(World &world, vector<Species> &species_list);
int collideMCC(Species &source, Field *target_den, double dt);
double getSigma(double v_mag);

/*random number generator*/
std::random_device rd;
std::mt19937 mt_gen(rd());		/*seed*/
std::uniform_real_distribution<double> rnd_dist(0, 1.0);
double rnd()
{
	return rnd_dist(mt_gen);
}

/*simulation parameters*/
const double PLASMA_DEN = 1e10;			//reference plasma density
const double kTe = 2;						//electron temperature;
const double DRIFT_VELOCITY = 7000;	
const double NEUTRAL_DEN = 1e13;
const int MCC_FREQ = -1;		/*how often to perform collisions*/
const int NUM_TS = 400;


/*globals*/
int mpi_rank=0;

/*helper function to print message only on root and abort*/
void error(string msg)
{
	if (mpi_rank==0) cerr<<msg<<endl;
	exit(-1);
}

/*main function*/
int main(int n_args, char *args[])
{
	int num_threads = 1;	//default value
	unsigned int max_threads = thread::hardware_concurrency();
	int doms[3]={1,1,1};

	 /* Initialize MPI */
	MPI_Init(NULL, NULL);

	/*get number of processes from MPI*/
	int mpi_size;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	/*get our rank*/
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	/*first argument is number of threads*/
	if (n_args>1)
	{
		num_threads = atoi(args[1]);
	}
	/*second argument is domain decomposition in format ni:nj:nk*/
	if (n_args>2)
	{
	    stringstream ss(args[2]);
	    vector<string> pieces;
	    string piece;
		while (getline(ss, piece, ':')) {
		        pieces.push_back(piece);
		    }
		if (pieces.size()!=3) error("Need to specify ni:nj:nk decomposition");
		/*conver to int*/
		doms[0] = stoi(pieces[0]);
		doms[1] = stoi(pieces[1]);
		doms[2] = stoi(pieces[2]);
		if (doms[0]*doms[1]*doms[2]!=mpi_size)
				error("Incorrect MPI size");
	}

	/*capture starting time*/
    auto clock_start = chrono::high_resolution_clock::now();

	/*initialize domain with global parameters*/
    World world(21,21,41);
    world.setSpacing(1e-2, 1e-2, 1e-2);
    world.setOrigin(-0.1, -0.1, 0);
    world.setDt(1e-7);

    /*this function will set our local domain based on the rank and domain decomposition*/
    world.initMPIDomain(mpi_rank,mpi_size,doms);

	/*set objects*/
    world.AddInlet();
    world.AddSphere(0,0,0.15,0.05,-100);

    if (mpi_rank==0)
    {
    	cout<<"Running with "<<num_threads<<" of "<<max_threads<<" maximum threads and "<<mpi_size<<" MPI domains"<<endl;
    }

    /*set number of threads*/
    world.setNumberOfThreads(num_threads);

	/*data structure to hold particle species*/
	vector<Species> species_list;

	/*add singly and doubly charged ions, and neutrals*/
	species_list.reserve(3);	/*avoid reallocation of memory and destructor calls*/
    species_list.emplace_back(world, "O+", 16*AMU, 1*QE, 5e1,(int)1e6);  /*singly charged ions*/
  //  species_list.emplace_back(world, "O++", 16*AMU, 2*QE, 5e1,(int)1e5); /*double charged ions*/
 //   species_list.emplace_back(world, "O", 16*AMU, 0, 1e5,(int)1e5);
	//Species *neutrals = &species_list[2];	        
	Species *neutrals =nullptr;

    /*add sources*/
	vector<Source> source_list;
	source_list.reserve(3);
    source_list.emplace_back(world,species_list[0],1*PLASMA_DEN,DRIFT_VELOCITY,0);
 //   source_list.emplace_back(world,species_list[1],0.1*PLASMA_DEN,DRIFT_VELOCITY,0);
  //  source_list.emplace_back(world,species_list[2],NEUTRAL_DEN,DRIFT_VELOCITY,0);

    /*solve initial potential*/
    PotentialSolver solver(world,SolverType::GS);
    solver.setElectronRef(0,PLASMA_DEN, kTe);
    solver.solvePotential();

    /*obtain electric field*/
    solver.computeEF();

	/*setup output handler*/
	Output output(world);
    
	/*flag whether we reached steady state*/
	bool steady_state = false;
	int np_sum_old = 0;

    /*show some diagnostics*/
    double lambda_d = sqrt(EPS_0*kTe/(PLASMA_DEN*QE));
    if (mpi_rank==0)
    	cout<<"Debye length = "<<lambda_d<<",\t Cell Spacing = "<<world.dh[0]<<endl;

    /* MAIN LOOP*/
    int part_tot = 0;
	int ts;
    for (ts = 1; ts<=NUM_TS; ts++)
    {
		/*inject particles*/
        for (Source &source:source_list)
            source.sample();

        for (Species &species:species_list) part_tot+=species.np;

        /*perform collisions, simulates collisions with neutrals*/
       // world.num_colls = 0;                
       // if (MCC_FREQ>0 && ts%MCC_FREQ==0)
        //    for (Species &species:species_list)
         //       world.num_colls = collideMCC(species,neutrals->den,world.dt*MCC_FREQ);

		/*move particles*/
		for (Species &species:species_list)
            species.move();

		/*compute charge density*/
        computeRhoi(world, species_list);

        /*update potential*/
        solver.solvePotential();

        /*obtain electric field*/
        solver.computeEF();


        /*every 10 time steps get global particle counts*/
        if (ts%10==0)
        {
        	int ns = species_list.size();
        	int *sp_count_local = new int[ns];
        	int *sp_count_global = new int[ns];

        	for (int i=0;i<ns;i++)
        		sp_count_local[i] = species_list[i].np;

        	/*sum up species counts over processors*/
        	MPI_Allreduce(sp_count_local,sp_count_global,ns,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

        	/*screen output*/
			if (mpi_rank==0)
			{
				cout<<ts<<":";
				/*output species info*/
				for (int i=0;i<ns;i++)
					cout<<"  "<<species_list[i].name<<":"<<sp_count_global[i];
//					cout<<"\t n_colls:"<<world.num_colls<<endl;
				cout<<endl;

				/*save 1D log data*/
				output.outputLog(ts, world.num_colls, species_list);
			}

	        /*check for steady state*/
	        if (!steady_state)
	        {
	        	int np_sum=0;
	        	for (int i=0;i<ns;i++)
	        		np_sum+=sp_count_global[i];

	            if (abs((np_sum-np_sum_old)/(double)np_sum)<0.005)
	            {
	                if(mpi_rank==0) cout<<"*** Reached Steady State ***"<<endl;
	                steady_state=true;
	            }

	            /*update old data*/
	            np_sum_old=np_sum;
	        }

			delete[] sp_count_local;
			delete[] sp_count_global;
        }

        /*start averaging at steady state*/
        if (steady_state && ts%10==0)
        {
            for (Species &species:species_list)
            {
                species.den_ave->updateAverage(species.den);
				for (int i=0;i<3;i++) species.vel_ave.f[i]->updateAverage(species.vel.f[i]);
            }
        }
            
        /*save results every 25 steps*/
        if (ts%25==0) 
        {			
            output.outputField(ts,species_list);
            output.outputTracer(ts,species_list[0],10);
            output.outputScatter(ts,species_list[0],10);
        }			

        world.time+=world.dt;
    }

    /*grab final time*/
 	auto clock_end = chrono::high_resolution_clock::now();
	std::chrono::duration<float> delta = clock_end-clock_start;

	/*save final results*/
	output.outputField(ts>0?ts-1:ts,species_list);

	if (mpi_rank==0)
    {
		output.outputPVD();
		cout << "Simulation took "<<delta.count()<< "s\n";
		cout<<"Done!"<<endl;
    }

	MPI_Finalize();

	return 0;
}

/*computes charge density of ions by multiplying number density by charge*/
void computeRhoi(World &world, vector<Species> &species_list)
{
    world.rhoi->clear();
	for (Species &species:species_list)
		for (int i=0;i<world.ni;i++)
			for (int j=0;j<world.nj;j++)
				for (int k=0;k<world.nk;k++)
					world.rhoi->data[i][j][k] += species.charge*species.den->data[i][j][k];
}	

/*demo of MCC collisions*/
int collideMCC(Species &species, Field *target_den, double dt)
{
	int num_colls = 0;
    for (int p=0;p<species.np;p++)
    {
    	Particle &part = species.particles[p];
        /*compute cross-section, function of relative velocity*/
        double v_mag = sqrt(part.vel[0]*part.vel[0]+
			            part.vel[1]*part.vel[1]+part.vel[2]*part.vel[2]);
            
        /*get total cross-section - normally you would add up different individual
        cross-sections            */
        double sigma = getSigma(v_mag);
            
        /*get target density*/
        double lc[3];
		target_den->world.XtoL(lc,part.pos);
        double nn = target_den->gather(lc);
            
        /*compute probability*/
        double P = 1 - exp(-nn*sigma*v_mag*dt);
            
        /*compare to a random number to see if collision happened*/
        if (P>=rnd())
        {
            /*for demo, approximate CEX*/
            /*need to sample a target particle*/
            double n[3];
            double theta = 2*PI*rnd();

            /*pick a random direction for n[2]*/
            double R = -1.0+2*rnd();
            double a = sqrt(1-R*R);

            n[0] = cos(theta)*a;
            n[1] = sin(theta)*a;
            n[2] = R;
                
            /*sample from Maxwellian*/
            double v_thermal = 300; /*should get from neutrals*/
            double v_target = v_thermal*2*(rnd()+rnd()+rnd()-1.5);
                
            for (int i=0;i<3;i++) part.vel[i] = v_target*n[i];
            num_colls++;
        }
    }
	return num_colls;
}
    
/*contant cross-section for the demo*/
double getSigma(double v_mag)
{
    return 1e-18;   /*in m^2*/
}
