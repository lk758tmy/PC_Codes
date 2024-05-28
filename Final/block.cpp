#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <math.h>
using namespace std;
using namespace std::chrono;
const int N=2500;
double min(double a,double b){
	return (a<b)?a:b;
}
int main(int argc,char *argv[]){
	int numProcess,myRank,npX,npY,lX,lY,myLX,myLY,src;
	int kStartX,kEndX,kStartY,kEndY,myIdX,myIdY;
	double *data=0,*result=0,*dataTmp=0,*resultTmp=0;
	int *sendSize1=0,*sendPos1=0,*sendSize2=0,*sendPos2=0;
	time_point<steady_clock> tStart,tEnd;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
	MPI_Comm_size(MPI_COMM_WORLD,&numProcess);
	
	if(numProcess==1) npX=npY=1;
	else if(numProcess==4) npX=npY=2;
	else if(numProcess==6) npX=2,npY=3;
	else if(numProcess==8) npX=2,npY=4;
	else if(numProcess==9) npX=3,npY=3;
	else{
		if(myRank==0) printf("Error: np should be 4/6/8/9.\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	/*npX=(int)sqrt(numProcess); npY=npX;
	if(npX*npY!=numProcess){
		if(myRank==0)
			printf("Error: np should be a square number\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}*/
	
	lX=N/npX; lY=N/npY;
	myIdX=myRank%npX; myIdY=myRank/npX; 
	kStartX=myIdX*lX; kStartY=myIdY*lY;
	if(myIdX+1==npX){
		kEndX=N; myLX=N-kStartX;
	}else{
		kEndX=(myIdX+1)*lX; myLX=lX;
	}
	if(myIdY+1==npY){
		kEndY=N; myLY=N-kStartY;
	}else{
		kEndY=(myIdY+1)*lY; myLY=lY;
	}

	MPI_Comm COMM_X,COMM_Y;
	MPI_Comm_split(MPI_COMM_WORLD,myIdX,myIdY,&COMM_X);
	MPI_Comm_split(MPI_COMM_WORLD,myIdY,myIdX,&COMM_Y);
	
	double *A=(double *)malloc(myLX*myLY*sizeof(double));
	double *D=(double *)malloc(myLX*myLY*sizeof(double));
	double *BX=(double *)malloc(myLY*sizeof(double));
	double *tmpY=(double *)malloc(myLX*sizeof(double)),*BY=0;
	
	//int myRX,myRY;
	//MPI_Comm_rank(COMM_X,&myRX); MPI_Comm_rank(COMM_Y,&myRY);
	//printf("%d %d %d %d %d %d %d %d %d %d %d\n",
	//	myRank,myIdX,myIdY,myLX,myLY,myRX,myRY,kStartX,kEndX,kStartY,kEndY);

	if(myRank==0){
		data=(double *)malloc(N*N*sizeof(double));
		result=(double *)malloc(N*N*sizeof(double));
		for(int i=0;i<N;i++)
			for(int j=0;j<N;j++)
				*(data+i*N+j)=(i%3+1)*(j%5+1)*5+((i+j)%7)+((j-i)/2)%3;
		
		tStart=steady_clock::now();
		
		sendPos1=(int *)malloc(npY*sizeof(int));
		sendSize1=(int *)malloc(npY*sizeof(int));
		for(int i=0;i<npY;i++){
			*(sendSize1+i)=lY*N; *(sendPos1+i)=i*lY*N;
		}
		*(sendSize1+npY-1)=N*N-*(sendPos1+npY-1);
	}
	if(myIdX==0){
		dataTmp=(double *)malloc(N*myLY*sizeof(double));
		resultTmp=(double *)malloc(N*myLY*sizeof(double));
		
		MPI_Scatterv(data,sendSize1,sendPos1,MPI_DOUBLE,dataTmp,N*myLY,MPI_DOUBLE,0,COMM_X);

		sendPos2=(int *)malloc(npX*sizeof(int));
		sendSize2=(int *)malloc(npX*sizeof(int));
		for(int i=0;i<npX;i++){
			*(sendSize2+i)=lX; *(sendPos2+i)=i*lX;
		}
		*(sendSize2+npX-1)=N-*(sendPos2+npX-1);	
	}
	for(int k=0;k<myLY;k++)
		MPI_Scatterv(dataTmp+k*N,sendSize2,sendPos2,MPI_DOUBLE,A+k*myLX,myLX,MPI_DOUBLE,0,COMM_Y);
	
	for(int i=0;i<myLY;i++)
		for(int j=0;j<myLX;j++)
			*(D+i*myLX+j)=*(A+i*myLX+j);
	
	for(int k=0;k<N;k++){
		
		if(k>=kStartY && k<kEndY){
			BY=D+myLX*(k-kStartY);
		}else BY=tmpY;
		src=min(k/lY,npY-1);
		MPI_Bcast(BY,myLX,MPI_DOUBLE,src,COMM_X);

		if(k>=kStartX && k<kEndX)
			for(int t=0;t<myLY;t++)
				*(BX+t)=*(D+(k-kStartX)+t*myLX);
		src=min(k/lX,npX-1);
		MPI_Bcast(BX,myLY,MPI_DOUBLE,src,COMM_Y);
		
		//MPI_Barrier(MPI_COMM_WORLD);
		
		for(int i=0;i<myLY;i++){
			if(i+kStartY==k) continue;
			for(int j=0;j<myLX;j++){
				if(j+kStartX==k) continue;
				if(*(D+i*myLX+j)>(*(BX+i)+*(BY+j)))
					*(D+i*myLX+j)=(*(BX+i)+*(BY+j));
			}
		}
	}
	
	for(int k=0;k<myLY;k++)
		MPI_Gatherv(D+k*myLX,myLX,MPI_DOUBLE,resultTmp+k*N,sendSize2,sendPos2,MPI_DOUBLE,0,COMM_Y);
	if(myIdX==0)
		MPI_Gatherv(resultTmp,myLY*N,MPI_DOUBLE,result,sendSize1,sendPos1,MPI_DOUBLE,0,COMM_X);
	
	if(myRank==0){
		tEnd=steady_clock::now();
		printf("%.3fs\n",(tEnd-tStart).count()/1000000000.0);
		printf("%d ",(int)*(result+472*N+1121));
		printf("%d ",(int)*(result+2253*N+768));
		printf("%d ",(int)*(result+1936*N+193));
		printf("%d ",(int)*(result+547*N+2038));
		printf("%d\n",(int)*(result+2353*N+1860));
		free(data); free(result); free(sendPos1); free(sendSize1);
	}
	if(myIdX==0){
		free(dataTmp); free(resultTmp); free(sendPos2); free(sendSize2);
	}
	free(A); free(D); free(tmpY); free(BX);
	
	MPI_Finalize();
	return 0;
}

//计算结果还不太对,且仅在npy!=npx时出现？
