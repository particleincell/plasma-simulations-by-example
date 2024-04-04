#include <mpi.h>
#include <iostream>
using namespace std;

int main(int n_args, char *args[]) 
{
 	MPI_Init(&n_args, &args);
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

 	int other_rank;
	if (my_rank == 0) {
		MPI_Ssend(&my_rank, 1, MPI_INT, 1, 42, MPI_COMM_WORLD);
		MPI_Recv(&other_rank, 1, MPI_INT, 1, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else if (my_rank == 1) {
		MPI_Ssend(&my_rank, 1, MPI_INT, 0, 42, MPI_COMM_WORLD);
		MPI_Recv(&other_rank, 1, MPI_INT, 0, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
	}
	MPI_Barrier(MPI_COMM_WORLD);

	cout<<"my_rank = "<<my_rank<<", other_rank = "<<other_rank<<endl;
	
	MPI_Finalize();
 	return 0;
}
