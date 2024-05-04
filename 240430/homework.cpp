#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi/mpi.h>
using namespace std;
const int N=1;
int main(int argc, char* argv[]){
	int myRank,numProcess,msgSize=sizeof(int)*N;
	
	MPI_Init(&argc,&argv);
	
	MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
	MPI_Comm_size(MPI_COMM_WORLD,&numProcess);
	int* sendMsg=(int *)malloc(msgSize);
	*sendMsg=myRank;
	int* recvMsg=(int *)malloc(msgSize*numProcess);
	
/*
Src	Dest
0:	1	2	3	4	5	6	7
1:	2	3	4	5	6	7	0	
2:	3	4	5	6	7	0	1
3:	4	5	6	7	0	1	2
4:	5	6	7	0	1	2	3
5:	6	7	0	1	2	3	4
6:	7	0	1	2	3	4	5
7:	0	1	2	3	5	6	7
In each layer, half of the processes will receive first
	and the other half will send first.
*/
	int dest=myRank,src=myRank,tag=0,tmp1,tmp2;
	MPI_Status stat;
	if(src<0) src=numProcess-1;
	if(dest>=numProcess) dest=0;
	
	*(recvMsg+myRank)=*sendMsg;
	for(int i=1;i<numProcess;i++){
		src--;	src=((src==-1)?(numProcess-1):src);
		dest++;	dest=((dest==numProcess)?0:dest);
		tmp1=i;	tmp2=myRank;
		while(tmp1%2==0){tmp1/=2; tmp2/=2; }
		if(tmp2%2==0){
			MPI_Send(sendMsg,N,MPI_INT,dest,tag,MPI_COMM_WORLD);
			//printf("Thread No.%d send to %d.\n",myRank,dest);
			MPI_Recv(recvMsg+src*N,N,MPI_INT,src,tag,MPI_COMM_WORLD,&stat);
		}else{
			MPI_Recv(recvMsg+src*N,N,MPI_INT,src,tag,MPI_COMM_WORLD,&stat);
			MPI_Send(sendMsg,N,MPI_INT,dest,tag,MPI_COMM_WORLD);
			//printf("Thread No.%d send to %d.\n",myRank,dest);
		}
		//sleep(1);
	}		

	printf("No.%d:\t",myRank);
	for(int i=0;i<numProcess;i++) printf("%d ",*(recvMsg+i));
	printf("\n");
	
	MPI_Finalize();
	return 0;
}
//mpirun -np 8 ./homework
