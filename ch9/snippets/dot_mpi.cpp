#include <mpi.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <sstream>

using namespace std;
using dvector = vector<double>;

int main(int n_args, char *args[]) {
  int mpi_size;
  int mpi_rank;

  MPI_Init(&n_args, &args);

  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  cout<<"I am "<<mpi_rank<<" of "<<mpi_size<<endl;
	
  auto start = chrono::system_clock::now();
  
  size_t n0 = 100000000;	//create two large arrays
  
  //get my limits
  size_t i1 = mpi_rank*n0/mpi_size;
  size_t i2 = mpi_rank<mpi_size-1 ? (mpi_rank+1)*n0/mpi_size  : n0;
  size_t my_n = i2-i1; 
  
  double *a = new double[my_n];
  double *b = new double[my_n];
  
  //set some values
  for (size_t i=0;i<my_n;i++) {a[i] = (i1+i)/(double)n0; b[i]=2*a[i];}
  
  double my_dot = 0;  
  
  for (size_t i=0;i<my_n;i++) my_dot+=a[i]*b[i];

  auto end = chrono::system_clock::now();
  
  //collect data 
  double dot;
  MPI_Allreduce(&my_dot, &dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //output timing info
  if (mpi_rank==0) {
  	std::chrono::duration<double,std::milli> elapsed = end - start;
  	cout<<"MPI time: "<<elapsed.count()<<" (ms), dot: "<<dot<<endl;
  }

  MPI_Finalize();
  delete[] a;
  delete[] b; 
  return 0;
}
