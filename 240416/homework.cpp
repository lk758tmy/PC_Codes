#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi/mpi.h>
using namespace std;
int main(int argc, char* argv[]){
	int myRank,numProcess;
	
	MPI_Init(&argc,&argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
	MPI_Comm_size(MPI_COMM_WORLD,&numProcess);
	int* sendMessage=(int *)malloc(sizeof(int)*1024);
	int* recvMessage=(int *)malloc(sizeof(int)*1024);
	int dest,src,tag=0,delta;
	MPI_Status stat;

/*
	0
	1	(0)						+1
	3	2	(1)	(0)				+2
	7	6	5	4	(3)...		+4
	15	14	13	12	11...		+8
*/
	dest=src=myRank; delta=1;
	if(myRank==0){ //root
		delta=1;
		printf("No.0 start to send.\n");
	}else{ //others
		do{
			delta<<=1;
		}while(delta<=myRank);
		src-=(delta>>1);
		MPI_Recv(recvMessage,1024,MPI_INT,src,tag,MPI_COMM_WORLD,&stat);
		sleep(1); //用于展示实验效果 ，实际不需要此指令
		printf("No.%d received.\n",myRank);	
	}
	while(1){
		dest=myRank+delta; delta*=2;
		if(dest>=numProcess) break;
		MPI_Send(sendMessage,1024,MPI_INT,dest,tag,MPI_COMM_WORLD);
		sleep(1); //用于展示实验效果 ，实际不需要此指令
	}
	
	MPI_Finalize();
	return 0;
}
//mpirun -np XX ./homework
