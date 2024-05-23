#include <mpi/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <chrono>
using namespace std;
using namespace std::chrono;
const int npX=4,npY=4,intervalsX=100,intervalsY=100;
const int itmax=4000;
const double tol=1e-6;
int main(int argc, char *argv[]){
	int numProcess,myRank,idX,idY,lX,lY;
	int startX,startY,startColour;
	double hx=1.0/intervalsX,hy=1.0/intervalsY,d,dx,dy;
	time_point<steady_clock> tStart,tEnd;
    int left,right,above,below;
    double rtmp,res,err;
	MPI_Datatype borderX,borderY;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcess);
	
    if (npX*npY!=numProcess){
        if(myRank==0)
           printf("Error: npx*npy not equal to np!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    idX=myRank%npX; idY=myRank/npX;
	
	lX=(intervalsX-1)/npX;
	lY=(intervalsY-1)/npY; //size of local grid
	startX=idX*lX;  startY=idY*lY;
	//global coordinate of local (0,0)
	if(idX==npX-1) lX=intervalsX-startX+1;
	else lX+=2;
	if(idY==npY-1) lY=intervalsY-startY+1;
	else lY+=2;
	startColour=(startX+startY)%2; //colour of local (0,0)
		
    double *u, *ut, *unew, *f, *utmp=0;
	int matrixSize=lX*lY;
    u    = (double*)calloc(matrixSize,sizeof(double));
    ut   = (double*)calloc(matrixSize,sizeof(double));
    unew = (double*)calloc(matrixSize,sizeof(double));
    f    = (double*)calloc(matrixSize,sizeof(double));
	//all initals are zero

    if(!u | !ut | !f | !unew){
        if(myRank==0)
            printf("Out of memory!\n");
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

	//true solution: ut(x,y) = - (x^2 + y^2)/4
    for(int i=0,xi,yj;i<lX;i++){
        for(int j=0;j<lY;j++){
            xi=hx*(startX+i);
            yj=hy*(startY+j);
            *(ut+i*lY+j)=-0.25*(xi*xi+yj*yj);
        }
    }
	
	//set right hand side
	d=0.5*hx*hx*hy*hy/(hx*hx+hy*hy);
	dx=d/(hx*hx);
	dy=d/(hy*hy);
    for(int i=0;i<lX;i++)
        for(int j=0;j<lY;j++)
            *(f+i*lY+j)=d;
 
	//boundary condition
    if(idX==0)
        for(int j=0;j<lY;j++)
            *(u+j)=*(unew+j)=*(ut+j);
    if(idX==npX-1)
        for(int j=0;j<lY;j++)
			*(u+(lX-1)*lY+j)=*(unew+(lX-1)*lY+j)=*(ut+(lX-1)*lY+j);
    if(idY==0)
        for(int i=0;i<lX;i++)
			*(u+i*lY)=*(unew+i*lY)=*(ut+i*lY);
    if(idY==npY-1)
        for(int i=0;i<lX;i++)
            *(u+(i+1)*lY-1)=*(unew+(i+1)*lY-1)=*(ut+(i+1)*lY-1);

    tStart=steady_clock::now();

/***************************************************
    Red-Black Gauss-Seidel iteration
****************************************************/
	
	//compute left, right, above and below processes 
    if (idX>0) left=idY*npX+idX-1;
    else left=MPI_PROC_NULL;
    if (idX<npX-1) right=idY*npX+idX+1;
    else right=MPI_PROC_NULL;
	if (idY>0) below=(idY-1)*npX+idX;
	else below=MPI_PROC_NULL;
	if (idY<npX-1) above=(idY+1)*npX+idX;
	else above=MPI_PROC_NULL;

	//create two new mpi datatype
    MPI_Type_vector(lX-2,1,lY,MPI_DOUBLE,&borderX);
    MPI_Type_contiguous(lY-2,MPI_DOUBLE,&borderY);
    MPI_Type_commit(&borderX);
    MPI_Type_commit(&borderY);

	//begin iteration loop
    for(int iter=0;iter<(itmax*2);iter++){
        if(iter%2==1) //Red
		for(int i=1;i<lX-1;i++)
        for(int j=1;j<lY-1;j++)
		if((startColour+i+j+iter/2)%2==0)
	        *(unew+i*lY+j)=*(f+i*lY+j)+dy*(*(u+i*lY+j+1)+*(u+i*lY+j-1))
								+dx*(*(u+(i+1)*lY+j)+*(u+(i-1)*lY+j));
		else *(unew+i*lY+j)=*(u+i*lY+j);
		
		if(iter%2==0) //Black
		for(int i=1;i<lX-1;i++)
		for(int j=1;j<lY-1;j++)
		if((startColour+i+j+iter/2)%2==1)
			*(unew+i*lY+j)=*(f+i*lY+j)+dy*(*(u+i*lY+j+1)+*(u+i*lY+j-1))
								+dx*(*(u+(i+1)*lY+j)+*(u+(i-1)*lY+j));
		else *(unew+i*lY+j)=*(u+i*lY+j);

		//only check residual once in every 25 loops, otherwise too slow!
		if(iter%50!=49) goto skipResidualCheck;
        res=0.0;
        for(int i=1;i<lX-1;i++)
        for(int j=1;j<lY-1;j++){
            rtmp=*(unew+i*lY+j)-*(u+i*lY+j);
			//unew≈(f+dx...+dy...) --> (unew-u)≈(f-u+dx...+dy...)=residual
			//Equals to the residual of last loop approximately.
            if(res<fabs(rtmp)) res=fabs(rtmp);
        }
        rtmp=res/d;

        MPI_Allreduce(&rtmp,&res,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        if(myRank==0)
			printf("iter=%d,residual≈%.4e\n",(iter+1)/2,res);  
		if(res<tol) break;
		
		skipResidualCheck:
		
		//update approximation *by swaping pointer*
		utmp=unew; unew=u; u=utmp;

		//update inner boundary values
		//u(1:lx-1,1) -> below, u(1:lx-1,0) <- below
        MPI_Sendrecv(u+lY+1,1,borderX,below,111,
                     u+lY,  1,borderX,below,111, 
                     MPI_COMM_WORLD,&status);
		//u(1:lx-1,ly-2) -> above, u(1:lx-1,ly-1) <- above
        MPI_Sendrecv(u+2*lY-2,1,borderX,above,111,
                     u+2*lY-1,1,borderX,above,111,
                     MPI_COMM_WORLD,&status);
		//u(1,1:ly-1) -> left, u(0,1:ly-1) <- left
        MPI_Sendrecv(u+lY+1,1,borderY,left,111,
                     u+1   ,1,borderY,left,111,
                     MPI_COMM_WORLD,&status);
		//u(lx-2,1:ly-1) -> right, u(lx-1,1:ly-1) <- right
        MPI_Sendrecv(u+lY*(lX-2)+1,1,borderY,right,111,
                     u+lY*(lX-1)+1,1,borderY,right,111,
                     MPI_COMM_WORLD,&status);	
    }

    tEnd=steady_clock::now();
    if (myRank==0)
        printf("Time=%.4es\n",(tEnd-tStart).count()/1000000000.0);

	//check solution
    err = 0.0;
    for (int i=1;i<lX-1;i++)
    for (int j=1;j<lY-1;j++){
         rtmp=fabs(*(ut+i*lY+j)-*(u+i*lY+j));
         if(err<rtmp) err=rtmp;
    }
	rtmp=err;

    MPI_Reduce(&rtmp,&err,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if(myRank==0) 
        printf("||u-ut||_inf = %.4e\n", err);

	MPI_Type_free(&borderX); MPI_Type_free(&borderY);
    free(u); free(ut); free(unew); free(f);
    MPI_Finalize();
	return 0;
}
