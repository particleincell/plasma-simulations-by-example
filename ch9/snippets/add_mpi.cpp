#include <mpi.h>
#include <iostream>
using namespace std;


int main(int n_args, char *args[]) 
{
 	MPI_Init(&n_args, &args);

	int mpi_size;
    int mpi_rank;

	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	cout<<"I am "<<mpi_rank<<" of "<<mpi_size<<endl;

	//root sends data to the worker and prints the result
	if (mpi_rank == 0) {
		double buf[2] = {1,2};

		//send an array of two doubles to rank 1 with tag set to 42
		MPI_Send(buf, 2, MPI_DOUBLE, 1, 42, MPI_COMM_WORLD);

		//receive a single double from rank 1 and store in c
		double c;
		MPI_Recv(&c, 1, MPI_DOUBLE, 1, 42, MPI_COMM_WORLD,MPI_STATUS_IGNORE);		
		
		//output
		cout<<"c = "<<c<<endl;
	}
	//worker receives data from root, adds it, and returns the result
	else if (mpi_rank == 1) {
		double buf[2];
		//receive an array of two doubles from rank 0 and store in buf
		MPI_Recv(buf, 2, MPI_DOUBLE, 0, 42, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		//perform operation
		double c = buf[0] + buf[1];

		//send the result to rank 0
		MPI_Send(&c, 1, MPI_DOUBLE, 0, 42, MPI_COMM_WORLD);
	}
	
	MPI_Finalize();
 	return 0;
}
