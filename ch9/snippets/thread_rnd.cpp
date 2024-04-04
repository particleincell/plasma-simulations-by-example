#include <thread>
#include <iostream>
#include <random>

using namespace std;
class Rnd {
public:
	Rnd(int seed): mt_gen{seed}, rnd_dist{0,1.0} {}
	double operator() () {return rnd_dist(mt_gen);}
protected:
	std::mt19937 mt_gen;	    //random number generator
	std::uniform_real_distribution<double> rnd_dist;  //uniform distribution
};
Rnd rnd_glob(0);	//global generator

class Worker {
public:
  Worker(int id):rnd{id} {}
  void operator() (int n_samples, bool use_global) { 
    double sum = 0;
	if (use_global)
		for (int i=0;i<n_samples;i++) sum += rnd_glob();
	else
		for (int i=0;i<n_samples;i++) sum += rnd();
    ave = sum/n_samples;
  }
  double ave;
  Rnd rnd;
};

int main(int n_args, char *args[]) {
	int num_threads = 12;
	bool use_global = false;
	if (n_args>1 && args[1][0]=='G') use_global = true;

	vector<Worker> workers;
	for (int i=0;i<num_threads;i++)
		workers.emplace_back(i);

	auto start = chrono::system_clock::now();
	vector<thread> threads;
	for (int i=0;i<num_threads;i++)
		threads.emplace_back(std::ref(workers[i]),100000,use_global);

	for (thread &t:threads)
		t.join();
	auto end = chrono::system_clock::now();

  	std::chrono::duration<double,std::milli> elapsed = end - start;
  	cout<<"Time with a "<<(use_global?"global":"local")<<" rnd(): "<<elapsed.count()<<" (ms)"<<endl;

	//output values
	/* for (Worker &w:workers)	cout<<w.ave<<" ";   cout<<endl;    */
	return 0;
}

