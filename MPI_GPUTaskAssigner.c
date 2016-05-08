#include <stdio.h>
#include <string.h>
#include "mpi.h"






int main(int argc,char *argv[])
{
    int numprocs, myid, dst; 
	int  root_process = 0;
    MPI_Status status;
    char message[100];
    MPI_Init(&argc, &argv); 
    MPI_Comm_rank(MPI_COMM_WORLD, &myid); 
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	if(myid==0)
    {
		sprintf(message, "Greetings from process %d!",myid);
		
		//read csv files => calculate nnz => sort by nnz => assign task to GPU devs		
		vector<>
		
		for (dst = 1; dst < numprocs; dst++) {
			MPI_Send(message,strlen(message)+1, MPI_CHAR, dst, 99, MPI_COMM_WORLD);
		} 
	}
	else 
    {   /* myid != 0 */
        MPI_Recv(message, 100, MPI_CHAR, root_process, 99, MPI_COMM_WORLD, &status);
        printf("%s myid %d \n", message, myid);      		
    }
    MPI_Finalize();
    return 0;
}
