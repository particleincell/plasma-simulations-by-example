#include <thread>
#include <iostream>
#include <vector>
#include <chrono>
#include <sstream>

using namespace std;
using dvector = vector<double>;

void ThreadWorker(double *a, double *b, size_t i1, size_t i2, double *res)
{
	*res = 0;
	for (size_t i=i1;i<i2;i++) 
		(*res)+=a[i]*b[i];
	//stringstream ss;
	//ss<<i1<<" "<<i2<<" "<<i2-i1<<" "<<(*res)<<endl;
	//cout<<ss.str();
}

int main(int n_args, char *args[]) {
  size_t n = 100000000;	//create two large arrays
  double *a = new double[n];
  double *b = new double[n];
  
  //set some values
  for (size_t i=0;i<n;i++) {a[i] = i/(double)n; b[i]=2*a[i];}
  
  //*** multithreaded version ***
  int num_threads = 4;
  if (n_args>1) num_threads = atoi(args[1]);
  cout<<"Running with "<<num_threads<<" threads"<<endl;
  auto start_threads = chrono::system_clock::now();
  
  vector<double> result(num_threads);

  vector<thread> threads;
  for (int i=0;i<num_threads;i++)
  {
     int i1 = i*n/num_threads;
     int i2 = i<num_threads-1 ? (i+1)*n/num_threads : n;
     
     threads.emplace_back(ThreadWorker,a,b,i1,i2,&result[i]);
  }
  //wait for threads to finish
  for (thread &t:threads) t.join();

  double dot_threads = 0;
  for (int i=0;i<num_threads;i++) {dot_threads += result[i];}
  auto end_threads = chrono::system_clock::now();
    
  //output timing info
  std::chrono::duration<double,std::milli> elapsed_threads = end_threads - start_threads;
  cout<<"Threads time: "<<elapsed_threads.count()<<" (ms), dot: "<<dot_threads<<endl;
 
  delete[] a;
  delete[] b; 
  return 0;
}
