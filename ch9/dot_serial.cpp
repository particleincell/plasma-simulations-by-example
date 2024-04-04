#include <thread>
#include <iostream>
#include <chrono>

using namespace std;

int main(int n_args, char *args[]) {
  size_t n = 100000000;	//create two large arrays
  double *a = new double[n];
  double *b = new double[n];
  
  //set some values
  for (size_t i=0;i<n;i++) {a[i] = i/(double)n; b[i]=2*a[i];}
  
  //*** serial version ***
  auto start_serial = chrono::system_clock::now();
  double dot_serial = 0;
  for (size_t i=0;i<n;i++) dot_serial += a[i]*b[i];
  auto end_serial = chrono::system_clock::now();
  
 //output timing info
  std::chrono::duration<double,std::milli> elapsed_serial = end_serial - start_serial;
  cout<<"Serial time: "<<elapsed_serial.count()<<" (ms), dot: "<<dot_serial<<endl;
 
  delete[] a;
  delete[] b; 
  return 0;
}
