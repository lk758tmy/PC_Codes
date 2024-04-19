#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi/mpi.h>
using namespace std;
const int N=1024*1024;
int main(int argc, char* argv[]){
	int myRank,numProcess;
	
	MPI_Init(&argc,&argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
	MPI_Comm_size(MPI_COMM_WORLD,&numProcess);
	int* sendMessage=(int *)malloc(sizeof(int)*N);
	int* recvMessage=(int *)malloc(sizeof(int)*N);
	int dest=myRank+1,src=myRank-1,tag=0;
	MPI_Status stat;
	if(src<0) src=numProcess-1;
	if(dest>=numProcess) dest=0;
	
	
	if(myRank%2==0){
		MPI_Send(sendMessage,N,MPI_INT,dest,tag,MPI_COMM_WORLD);
		printf("Thread No.%d sent.\n",myRank);
		sleep(1);
		MPI_Recv(recvMessage,N,MPI_INT,src,tag,MPI_COMM_WORLD,&stat);
		printf("Thread No.%d received.\n",myRank);
	}else{
		MPI_Recv(recvMessage,N,MPI_INT,src,tag,MPI_COMM_WORLD,&stat);
		printf("Thread No.%d received.\n",myRank);
		sleep(1);
		MPI_Send(sendMessage,N,MPI_INT,dest,tag,MPI_COMM_WORLD);
		printf("Thread No.%d sent.\n",myRank);
	}
	
	MPI_Finalize();
	return 0;
}
//mpirun -np 8 ./homework
