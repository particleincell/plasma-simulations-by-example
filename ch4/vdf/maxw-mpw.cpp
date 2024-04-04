//demo of loading VDF from given function with variable mpw


#include <vector>
#include <chrono>
#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

using namespace std;

using dvector = vector<double>;

/*object for sampling random numbers*/
class Rnd {
public:
	//constructor: set initial random seed and distribution limits
	Rnd(): mt_gen{std::random_device()()}, rnd_dist{0,1.0} {}
	double operator() () {return rnd_dist(mt_gen);}

protected:
	std::mt19937 mt_gen;	    //random number generator
	std::uniform_real_distribution<double> rnd_dist;  //uniform distribution
};

Rnd rnd;

int main(int num_args, char* args[]) 
{
	const int NUM_BINS = 41;	//number of uniques we would like
	double PARTS_PER_BIN = 10;			//number of particles to load per bin
	double MPW0 = 1;					//nominal particle weight
	double np_uniform = 200;   //number of particles to load with const weight     
	double vth = 1e5;			//thermal velocity
	const double pi = acos(-1.0);

	dvector bins(NUM_BINS);

	double v_min = -6*vth;
	double v_max = 6*vth;
	double dv = (v_max-v_min)/(NUM_BINS);

	chrono::time_point<chrono::high_resolution_clock> time_start = chrono::high_resolution_clock::now();	//save starting time point
	
	double I = 0;
	double a = 1/(sqrt(pi)*vth);
	
	//loop over bins
	for (int i=0;i<NUM_BINS;i++)
	{
		//evaluate normalized distribution function
		double v = v_min + (i)*dv;
		double fm = a*exp(-v*v/(vth*vth));
		I += fm*dv;	
		//particle weight
		double mpw = (MPW0*np_uniform)*fm*dv/PARTS_PER_BIN;
		
	 	for (int p=0;p<PARTS_PER_BIN;p++) 
		{
		     //pick random velocity in this bin
		     v = v_min+i*dv + rnd()*dv;
		     
		     //bin result, add to nearest
		     int bin = (int)((v-v_min)/dv);
		     bins[bin] += mpw;
		}
	}

	auto time_now = chrono::high_resolution_clock::now();
  	chrono::duration<double,  std::nano> time_delta = time_now-time_start;
	cout<<"Time per sample: "<<time_delta.count()/(NUM_BINS*PARTS_PER_BIN)<<" ns"<<endl;

	//normalize bins
	double sp_sum = 0;
	for (int i=0;i<NUM_BINS;i++) 
		sp_sum += bins[i];

	cout<<"spwt_sum = "<<sp_sum<<endl;
	cout<<"I = "<<I<<endl;

	double max_val = 0;
	for (int i=0;i<NUM_BINS;i++) 
		if (bins[i]>max_val) max_val=bins[i];
	for (int i=0;i<NUM_BINS;i++) 
		bins[i]/=max_val;
		
	bins[0]*=2;
	bins[NUM_BINS-1]*=2;

	//write to a file
	ofstream out("maxw-mpw.csv");
	out<<"vel/vth,f_num,f_th\n";
	for (int i=0;i<NUM_BINS;i++) 
	{
		double v = v_min+i*dv;
		double f_th = exp(-v*v/(vth*vth));
		out<<v/vth<<","<<bins[i]<<","<<f_th<<"\n";
	}


}
