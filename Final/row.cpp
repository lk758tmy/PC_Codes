#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
using namespace std;
using namespace std::chrono;
const int N=2500;
double min(double a,double b){
	return (a<b)?a:b;
}
int main(int argc,char *argv[]){
	int numProcess,myRank,kStart,kEnd,l,myL,src;
	double *data=0,*result=0; int *sendSize=0,*sendPos=0;
	time_point<steady_clock> tStart,tEnd;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
	MPI_Comm_size(MPI_COMM_WORLD,&numProcess);

	l=N/numProcess; kStart=myRank*l;
	if(myRank+1==numProcess){
		kEnd=N; myL=N-kStart;
	}else{
		kEnd=(myRank+1)*l; myL=l;
	}
	
	double *A=(double *)malloc(myL*N*sizeof(double));
	double *D=(double *)malloc(myL*N*sizeof(double));
	double *tmp=(double *)malloc(N*sizeof(double)),*B=0;
	
	if(myRank==0){
		data=(double *)malloc(N*N*sizeof(double));
		result=(double *)malloc(N*N*sizeof(double));
		for(int i=0;i<N;i++)
			for(int j=0;j<N;j++)
				*(data+i*N+j)=(i%3+1)*(j%5+1)*5+((i+j)%7)+((j-i)/2)%3;
		
		tStart=steady_clock::now();
		
		sendPos=(int *)malloc(numProcess*sizeof(int));
		sendSize=(int *)malloc(numProcess*sizeof(int));
		for(int i=0;i<numProcess;i++){
			*(sendSize+i)=l*N; *(sendPos+i)=i*l*N;
		}
		*(sendSize+numProcess-1)=N*N-*(sendPos+numProcess-1);
	}
	MPI_Scatterv(data,sendSize,sendPos,MPI_DOUBLE,A,N*myL,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	for(int i=0;i<myL;i++)
		for(int j=0;j<N;j++)
			*(D+i*N+j)=*(A+i*N+j);
	
	for(int k=0;k<N;k++){
		if(k>=kStart && k<kEnd)
			B=D+N*(k-kStart);
		else B=tmp;
		src=min(k/l,numProcess-1);
		MPI_Bcast(B,N,MPI_DOUBLE,src,MPI_COMM_WORLD);
		
		for(int i=0;i<myL;i++){
			if(i+kStart==k) continue;
			for(int j=0;j<N;j++){
				if(j==k) continue;
				if(*(D+i*N+j)>(*(D+i*N+k)+*(B+j)))
					*(D+i*N+j)=(*(D+i*N+k)+*(B+j));
			}
		}
		//MPI_Barrier(MPI_COMM_WORLD);
	}
	
	MPI_Gatherv(D,myL*N,MPI_DOUBLE,result,sendSize,sendPos,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	if(myRank==0){
		tEnd=steady_clock::now();
		printf("%.3fs\n",(tEnd-tStart).count()/1000000000.0);
		printf("%d ",(int)*(result+472*N+1121));
		printf("%d ",(int)*(result+2253*N+768));
		printf("%d ",(int)*(result+1936*N+193));
		printf("%d ",(int)*(result+547*N+2038));
		printf("%d\n",(int)*(result+2353*N+1860));
		free(data); free(result); free(sendPos); free(sendSize);
	}
	free(A); free(D); free(tmp);
	
	MPI_Finalize();
	return 0;
}
