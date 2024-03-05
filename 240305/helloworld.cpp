#include <cstdio>
#include <mpi/mpi.h>
using namespace std;
int main(int argc, char *argv[]){
	int myrank;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank); //进程号
	printf("%d helloworld\n",myrank);
	MPI_Finalize();
	return 0;
}
// mpirun -np x ./helloworld    x为进程数
