#include <thread>
#include <iostream>
#include <chrono>

using namespace std;

double dot_serial(double *a, double *b, size_t n) {
  double dot = 0;
  for (size_t i=0;i<n;i++) dot += a[i]*b[i];
  return dot;
}

int main(int n_args, char *args[]) {
  size_t n = 100000000;	//create two large arrays
  double *a = new double[n];
  double *b = new double[n];
  
  //set some values
  for (size_t i=0;i<n;i++) {a[i] = i/(double)n; b[i]=2*a[i];}
  
  //*** serial version ***
  auto start = chrono::system_clock::now();
  double dot = dot_serial(a,b,n);
  auto end = chrono::system_clock::now(); 
  std::chrono::duration<double,std::milli> elapsed = end - start;
  cout<<"Serial time: "<<elapsed.count()<<" (ms), dot: "<<dot<<endl;
 
  delete[] a;
  delete[] b; 
  return 0;
}

