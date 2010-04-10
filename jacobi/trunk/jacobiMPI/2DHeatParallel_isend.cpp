// mpiheatP.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"


void heatParallel(int n, double epsilon)
{
	//for MPI
	int myid,num_processor;
	int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
	double startwtime = 0.0, endwtime;
	MPI_Status statusP,statusN;
	MPI_Request requestSP,requestSN,requestRP,requestRN;
	//for Compute
	int				i,j;
	int				goal=(n-2)*(n-2);
	int				num,totalnum;
	int				jobstart,jobend,task,jobleft;
	double			mean,averagemean;//Average boundary value
	double**		u;//Old values
	double**		w;//New values
	double**		temp;
	u= (double **)malloc(sizeof(double *)*n);
	for (int i=0; i<n; i++)
		u[i] = (double *)malloc(sizeof(double)*n);
	w= (double **)malloc(sizeof(double *)*n);
	for (int i=0; i<n; i++)
		w[i] = (double *)malloc(sizeof(double)*n);
	int step=0;
	//MPI Initialize
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&num_processor);
	MPI_Get_processor_name(processor_name,&namelen);
	//Set boundary values and compute mean boundary value
	task=(n-2)/num_processor;
	jobleft=n-2-num_processor*task;
	if (myid<jobleft)
	{
		jobstart=myid*(task+1)+1;
		jobend=(myid+1)*(task+1)+1;
	}
	else
	{
		jobstart=jobleft+myid*task+1;
		jobend=jobleft+(myid+1)*task+1;
	}	
	//Data Initialize
	mean=0.0;	
	for(i=1;i<n-1;i++)//由于四个角上的数值用不到，所以忽略四个角。
	{
		w[0][i]=u[0][i]=100.0;
		w[n-1][i]=u[n-1][i]=0.0;
		w[i][0]=w[i][n-1]=u[i][0]=u[i][n-1]=100.0;
		mean=mean+u[i][0]+u[i][n-1]+u[0][i]+u[n-1][i];	
	}
	averagemean=mean/(4.0*n-8);
	//Initialize interior values
	for(i=jobstart;i<jobend;i++)
		for(j=1;j<n-1;j++)
			u[i][j]=averagemean;
	if (myid == 0) 
	{
		fprintf(stdout, "N=%d  Epsilon=%f  Computing...  ",n,epsilon);
		fflush(stdout);
	}
	int done=0;
	startwtime=MPI_Wtime();	
	//Compute steady-state solution
	while (!done)
	{
		//Synchronize boundary data
		if (myid%2==0)
		{
			//Next
			if (myid!=num_processor-1)
			{
				MPI_Isend(u[jobend-1],n,MPI_DOUBLE,myid+1,0,MPI_COMM_WORLD,&requestSN);
				MPI_Irecv(u[jobend],n,MPI_DOUBLE,myid+1,1,MPI_COMM_WORLD,&requestRN);
			}
			//Previous
			if (myid!=0)
			{
				MPI_Isend(u[jobstart],n,MPI_DOUBLE,myid-1,1,MPI_COMM_WORLD,&requestSP);
				MPI_Irecv(u[jobstart-1],n,MPI_DOUBLE,myid-1,0,MPI_COMM_WORLD,&requestRP);				
			}
		}
		else
		{
			//Previous
			MPI_Irecv(u[jobstart-1],n,MPI_DOUBLE,myid-1,0,MPI_COMM_WORLD,&requestRP);
			MPI_Isend(u[jobstart],n,MPI_DOUBLE,myid-1,1,MPI_COMM_WORLD,&requestSP);
			//Next			
			if (myid!=num_processor-1)
			{				
				MPI_Irecv(u[jobend],n,MPI_DOUBLE,myid+1,1,MPI_COMM_WORLD,&requestRN);
				MPI_Isend(u[jobend-1],n,MPI_DOUBLE,myid+1,0,MPI_COMM_WORLD,&requestSN);
			}
		}		
		//Compute Inner Data 
		num=0;
		for(i=jobstart+1;i<jobend-1;i++)
		{
			for(j=1;j<n-1;j++)
			{
				w[i][j]=(u[i-1][j]+u[i+1][j]+u[i][j-1]+u[i][j+1])/4.0;
				if (fabs(w[i][j]-u[i][j])<=epsilon) 
					num++;
			}
		}
		//Compute boundary data
		if (myid==0)
		{
			i=jobstart;
			for(j=1;j<n-1;j++)
			{
				w[i][j]=(u[i-1][j]+u[i+1][j]+u[i][j-1]+u[i][j+1])/4.0;
				if (fabs(w[i][j]-u[i][j])<=epsilon) 
					num++;
			}
			MPI_Wait(&requestRN,&statusN);
			i=jobend-1;
			for(j=1;j<n-1;j++)
			{
				w[i][j]=(u[i-1][j]+u[i+1][j]+u[i][j-1]+u[i][j+1])/4.0;
				if (fabs(w[i][j]-u[i][j])<=epsilon) 
					num++;
			}
			MPI_Wait(&requestSN,&statusN);
		}
		else if (myid==num_processor-1)
		{
			i=jobend-1;
			for(j=1;j<n-1;j++)
			{
				w[i][j]=(u[i-1][j]+u[i+1][j]+u[i][j-1]+u[i][j+1])/4.0;
				if (fabs(w[i][j]-u[i][j])<=epsilon) 
					num++;
			}
			MPI_Wait(&requestRP,&statusP);
			i=jobstart;
			for(j=1;j<n-1;j++)
			{
				w[i][j]=(u[i-1][j]+u[i+1][j]+u[i][j-1]+u[i][j+1])/4.0;
				if (fabs(w[i][j]-u[i][j])<=epsilon) 
					num++;
			}
			MPI_Wait(&requestSP,&statusP);
		}
		else
		{
			//还可以MPI_Waitany()来提高性能***
			MPI_Wait(&requestRP,&statusP);
			i=jobstart;
			for(j=1;j<n-1;j++)
			{
				w[i][j]=(u[i-1][j]+u[i+1][j]+u[i][j-1]+u[i][j+1])/4.0;
				if (fabs(w[i][j]-u[i][j])<=epsilon) 
					num++;
				
			}
			MPI_Wait(&requestRN,&statusN);
			i=jobend-1;
			for(j=1;j<n-1;j++)
			{
				w[i][j]=(u[i-1][j]+u[i+1][j]+u[i][j-1]+u[i][j+1])/4.0;
				if (fabs(w[i][j]-u[i][j])<=epsilon) 
					num++;
			}			
			MPI_Wait(&requestSP,&statusP);
			MPI_Wait(&requestSN,&statusN);
		}
		totalnum=0;
		MPI_Allreduce(&num,&totalnum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
		if (totalnum==goal)
			done=1;
		step++;		
		temp = u;
		u = w;
		w = temp;		
	}	
	endwtime=MPI_Wtime();
	if (myid == 0) 
	{
		fprintf(stdout, "Done Computing.\n");
		fprintf(stdout, "N=%d  Epsilon=%f  Step=%d  Compute Time=%f\n\n",n,epsilon,step,endwtime-startwtime);
		fflush(stdout);		
	}
}

int main(int argc, char* argv[])
{
	
	MPI_Init(&argc,&argv);
	int				n;
	double			epsilon=0.01;

	n = 100;		heatParallel(n,epsilon);
	n = 500;		heatParallel(n,epsilon);
	n = 1000;		heatParallel(n,epsilon);
	n = 2000;		heatParallel(n,epsilon);
	n = 4000;		heatParallel(n,epsilon);
	n = 6000;		heatParallel(n,epsilon);
	n = 8000;		heatParallel(n,epsilon);

	MPI_Finalize();

	getchar();
	return 0;
}

   