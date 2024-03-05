#include <cstdio>
#include <cmath>
#include <mpi/mpi.h>
using namespace std;
int main(int argc, char* argv[]){
	int myRank,n,numProcess; double result=0;
	
	n=atoi(argv[1]);
	//scanf("%d",&n); 有问题但不知到为什么，据说进程再Init前已经创建了？
	MPI_Init(&argc,&argv);		

	MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
	MPI_Comm_size(MPI_COMM_WORLD,&numProcess); //进程个数
	double h=1.0/n,x,myResult=0;
	for(int i=myRank;i<n;i+=numProcess){
		x=(i+0.5)*h;
		myResult+=(4.0/(1+x*x));
	} // PI=4*integrate(0,1,1/(1+x*x));
	myResult*=h;
	//MPI_Bcast(&myResult,1,MPI_DOUBLE,0,MPI_COMM_WORLD); //广播
	MPI_Reduce(&myResult,&result,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); //归约
	  
	if(myRank==0) printf("A.Error:%.6e\n",abs(result-M_PI));
	
	MPI_Finalize();
	return 0;
}
//mpirun -np X ./pi n   X为进程数量 n为分段数量
