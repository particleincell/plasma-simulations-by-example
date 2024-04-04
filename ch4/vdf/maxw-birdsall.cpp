//demo of loading maxwellian VDF

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
	const int NUMS = 10000000;		//number of samples
	const int NUM_BINS = 2000;		//number of uniques we would like
	int M = 3;   		          	//parameter for Birdsall's method

	if (num_args>1) M = atoi(args[1]);
	if (M<1) M=1; else if (M>12) M=12;
	cout<<"Running with M = "<<M<<endl;

	double vth = 1e5;				//thermal velocity

	dvector bins(NUM_BINS);

	double v_min = -6*vth;
	double v_max = 6*vth;
	double dv = (v_max-v_min)/(NUM_BINS);

	chrono::time_point<chrono::high_resolution_clock> time_start = chrono::high_resolution_clock::now();	//save starting time point
	
	for (int s=0;s<NUMS;s++)
	{
		//sample Maxwellian using Birdsall's method
		double Rsum = 0;
		for (int i=0;i<M;i++)
			Rsum += rnd();
		
		double v = sqrt(0.5)*vth*(Rsum-M/2.0)/sqrt(M/12.0);
		
		//bin result, add to nearest
		int bin = (int)((v-v_min)/dv+0.5);
		bins[bin]++;
	}

	auto time_now = chrono::high_resolution_clock::now();
  	chrono::duration<double,  std::nano> time_delta = time_now-time_start;
	cout<<"Time per sample: "<<time_delta.count()/NUMS<<" ns"<<endl;

	//normalize bins
	double max_val = 0;
	for (int i=0;i<NUM_BINS;i++) 
		if (bins[i]>max_val) max_val=bins[i];

	for (int i=0;i<NUM_BINS;i++) 
		bins[i]/=max_val;
		
	//write to a file
	ofstream out("maxw-birdsall"+to_string(M)+".csv");
	out<<"vel/vth,f_num,f_th\n";
	for (int i=0;i<NUM_BINS;i++) 
	{
		double v = v_min+i*dv;
		double f_th = exp(-v*v/(vth*vth));
		out<<v/vth<<","<<bins[i]<<","<<f_th<<"\n";

	}
	return 0;
}
