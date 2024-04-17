#include <stdio.h>
#include <stdlib.h>
#include <chrono>
//使用C++中的chrono进行计时，其余均为C语言中的指令
#include <random>
#include <time.h>
#include <memory.h>
using namespace std;
using namespace std::chrono;
void rand_matrix(double *m,int n){
	srand(time(0));
	for(int i=0;i<n*n;i++)
		*(m+i)=double(rand()+1)/RAND_MAX;
	return ;
}
void dgemm_ijk(double *A,double *B,double *C,int n){
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			for(int k=0;k<n;k++)
				*(C+i*n+j)+=(*(A+i*n+k))*(*(B+k*n+j));
	return ;
}
void dgemm_jki(double *A,double *B,double *C,int n){
	for(int j=0;j<n;j++)
		for(int k=0;k<n;k++)
			for(int i=0;i<n;i++)
				*(C+i*n+j)+=(*(A+i*n+k))*(*(B+k*n+j));
	return ;
}
void dgemm_ikj(double *A,double *B,double *C,int n){
	for(int i=0;i<n;i++)
		for(int k=0;k<n;k++)
			for(int j=0;j<n;j++)
				*(C+i*n+j)+=(*(A+i*n+k))*(*(B+k*n+j));
	return ;
}
int main(){
	int matrix_size=1024*1024*sizeof(double),N=20;
	time_point<steady_clock> t0,t1,t2,t3;
	double *A=(double *)malloc(matrix_size);
	double *B=(double *)malloc(matrix_size);
	double *C=(double *)malloc(matrix_size);
	double t_ijk=0,t_jki=0,t_ikj=0;
	
	for(int i=0;i<N;i++){
		printf("Test %d\n",i+1);
		rand_matrix(A,1024);
		rand_matrix(B,1024);

		t0=steady_clock::now();
		memset(C,0,matrix_size);
		dgemm_ijk(A,B,C,1024);
		t1=steady_clock::now();
		memset(C,0,matrix_size);
		dgemm_jki(A,B,C,1024);
		t2=steady_clock::now();
		memset(C,0,matrix_size);
		dgemm_ikj(A,B,C,1024);
		t3=steady_clock::now();
		
		t_ijk+=((t1-t0).count()/1000000000.0);
		t_jki+=((t2-t1).count()/1000000000.0);
		t_ikj+=((t3-t2).count()/1000000000.0);
	}
	
	printf("\tTime(s)\tAvg.Time(s)\n");
	printf("i-j-k\t%.3f\t%.3f\n",t_ijk,t_ijk/N);
	printf("j-k-i\t%.3f\t%.3f\n",t_jki,t_jki/N);
	printf("i-k-j\t%.3f\t%.3f\n",t_ikj,t_ikj/N);
	free(A); free(B); free(C);
	return 0;
}
