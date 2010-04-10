
#include "stdafx.h"

//*****************************************************************************
void getLocalJob(const int n, 
				 const int nodeCount, const int nodeIndex, 
				 int *jobStartingPoint, int *jobEnd)
{
	int jobSize = (n - 2) / nodeCount;
	int jobLeft = (n - 2) % nodeCount;
	if (nodeIndex < jobLeft)
	{
		*jobStartingPoint	= nodeIndex * (jobSize + 1) + 1;
		*jobEnd				= (myid + 1) * (task + 1) + 1;
	}
	else
	{
		*jobStartingPoint	= jobLeft + nodeIndex * jobSize + 1;
		*jobEnd				= jobLeft + (nodeIndex + 1) * jobSize + 1;
	}
	return;
}

//*****************************************************************************
void matrix1DInit(const int n, const struct boundary b, 
				  const int jobStartingPoint, const int jobEnd,
				  double *m, double *w)
{
	for(int i = 0; i < n - 1; i++)
	{
		m[i * n + n]			= w[i * n + n]			= b.left; 
		m[i * n + n - 1]		= w[i * n + n - 1]		= b.right; 
		m[i]					= w[i]					= b.up; 
		m[n * n - n + i + 1]	= w[n * n - n + i + 1]	= b.down;
	}
	for(int i = jobStartingPoint; i < jobEnd; i++)
		for(int j = 1; j < n - 1; j++)
			m[i * n + j] = w[i * n + j] = b.averageValue;
	return;
}

void rowDataIteration(const int n, const long step
					  const double epsilon, const int rowIndex,
					  double *m, double *w, int *localCount)
{
	for(int j = 1; j < n - 1; j++)
	{
		w[rowIndex * n + j] = 
			(m[rowIndex * n + j - n] + m[rowIndex * n + j + n]
			 + m[rowIndex * n + j - 1] + m[rowIndex * n + j + 1]) / 4.0;
		if (step % JUMP == 0)
			if(fabs(w[rowIndex * n + j] - m[rowIndex * n + j]) < epsilon) 
				(*localCount)++;
	}
}

//*****************************************************************************
void jacobiMPIIterationEpsilon_1D(const int n, const double epsilon, 
								  long *step, const struct boundary b, 
								  double *m, double *w,
								  double *initTime, double *iterTime)
{
	//MPI basic parameters
	int				nodeIndex;
	int				nodeCount;
	int				nodeNameLength;
    char			nodeName[MPI_MAX_PROCESSOR_NAME];
	//MPI Initialize
	MPI_Comm_rank(MPI_COMM_WORLD, &nodeIndex);
	MPI_Comm_size(MPI_COMM_WORLD, &nodeCount);
	MPI_Get_processor_name(nodeName, &nodeNameLength);

	double			*temp;
	
	//timer
	LARGE_INTEGER	nStartCounter, nStopCounter;

	//init data
	printf("--Data initing(%d, %lf).....", n, epsilon);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	//set job for local node
	int				jobStartingPoint, jobEnd;
	getLocalJob(n, nodeCount, nodeIndex, &jobStartingPoint, &jobEnd);
	matrix1DInit(n, b, jobStartingPoint, jobEnd, m, w); 
	//timer ends
	QueryPerformanceCounter(&nStopCounter);
	*initTime = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");
	//iteration
	printf("--Computing(%d, %lf).....", n, epsilon);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	//MPI Synchronous parameters
	MPI_Status		statusPrevious, statusNext;
	MPI_Request		requestSendPrevious, requestSendNext,
	MPI_Request		requestReceivePrevious, requestReceiveNext;
	//iteration parameters
	int				localCount, totalCount;
	int				countGoal = (n - 2) * (n - 2);
	*step = 0;
	while(totalCount < countGoal)
	{
		(*step)++;

		//Synchronize boundary data
		if (nodeIndex % 2 == 0)
		{
			//Next
			if (nodeIndex != nodeCount - 1)
			{
				MPI_Isend(m[jobEnd - 1], n, MPI_DOUBLE, 
					nodeIndex + 1, 0, MPI_COMM_WORLD, &requestSendNext);
				MPI_Irecv(m[jobEnd], n, MPI_DOUBLE, 
					nodeIndex + 1, 1, MPI_COMM_WORLD, &requestReceiveNext);
			}
			//Previous
			if (nodeIndex != 0)
			{
				MPI_Isend(m[jobStartingPoint], n, MPI_DOUBLE,
					nodeIndex - 1, 1, MPI_COMM_WORLD, &requestSendPrevious);
				MPI_Irecv(m[jobStartingPoint - 1], n , MPI_DOUBLE,
					nodeIndex - 1, 0, MPI_COMM_WORLD, &requestReceivePrevious);				
			}
		}
		else
		{
			//Previous
			MPI_Irecv(m[jobStartingPoint - 1], n, MPI_DOUBLE, 
				nodeIndex - 1, 0, MPI_COMM_WORLD, &requestReceivePrevious);
			MPI_Isend(m[jobStartingPoint], n, MPI_DOUBLE, 
				nodeIndex - 1, 1, MPI_COMM_WORLD, &requestSendPrevious);
			//Next			
			if (nodeIndex != nodeCount - 1)
			{				
				MPI_Irecv(m[jobEnd], n, MPI_DOUBLE, 
					nodeIndex + 1, 1, MPI_COMM_WORLD, &requestSendNext);
				MPI_Isend(m[jobEnd - 1], n, MPI_DOUBLE, 
					nodeIndex + 1, 0, MPI_COMM_WORLD, &requestReceiveNext);
			}
		}
		//Compute Inner Data 
		localCount = 0;		
		for(int i = jobStartingPoint + 1; i < jobEnd - 1; i++)
			rowDataIteration(n, *step, epsilon, i, m, w, &localCount);
		//Compute boundary data
		int rowIndex;
		if (nodeIndex == 0)
		{
			rowIndex = jobStartingPoint;
			rowDataIteration(n, *step, epsilon, rowIndex, m, w, &localCount);
			MPI_Wait(&requestReceiveNext, &statusNext);
			rowIndex = jobEnd - 1;
			rowDataIteration(n, *step, epsilon, rowIndex, m, w, &localCount);
			MPI_Wait(&requestSendNext, &statusNext);
		}
		else if (myid==num_processor-1)
		{
			rowIndex = jobEnd - 1;
			rowDataIteration(n, *step, epsilon, rowIndex, m, w, &localCount);
			MPI_Wait(&requestReceivePrevious, &statusPrevious);
			rowIndex = jobStartingPoint;
			rowDataIteration(n, *step, epsilon, rowIndex, m, w, &localCount);
			MPI_Wait(&requestSendPrevious, &statusPrevious);
		}
		else
		{
			//还可以MPI_Waitany()来提高性能***
			MPI_Wait(&requestReceivePrevious, &statusPrevious);
			rowIndex = jobStartingPoint;
			rowDataIteration(n, *step, epsilon, rowIndex, m, w, &localCount);
			MPI_Wait(&requestReceiveNext, &statusNext);
			rowIndex = jobEnd - 1;
			rowDataIteration(n, *step, epsilon, rowIndex, m, w, &localCount);			
			MPI_Wait(&requestSendPrevious, &statusPrevious);
			MPI_Wait(&requestSendNext, &statusNext);
		}
		totalCount = 0;
		MPI_Allreduce(&localCount, 
			&totalCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		temp = m; m = w; w = temp;
	}
	//timer ends
	QueryPerformanceCounter(&nStopCounter);
	*iterTime = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");

	return;
}

//*****************************************************************************
void jacobiMPI_1D(int n, double epsilon, 
				  long step, struct boundary b, char *outFile)
{
	printf("Jacobi MPI 1D -\n");
	printf("--n=%d, e=%lf, step=%ld\n--LURD: %lf, %lf, %lf, %lf\n",
		n, epsilon, step, b.left, b.up, b.right, b.down);
	//more paramenters
	double			*m = (double *)malloc(sizeof(double) * n * n);
	double			*w = (double *)malloc(sizeof(double) * n * n);

	//timer
	LARGE_INTEGER	nStartCounter, nStopCounter;
	double			nTime1, nTime2, nTime3;

	//jacobi serial 1D solution
	if (epsilon != 0)
	{
		printf("--Epsilon mode\n");
		jacobiMPIIterationEpsilon_1D(n, epsilon, &step, 
										b, m, w, &nTime1, &nTime2);
		printf("--Step = %ld\n", step);
	}
	else 
	{
		printf("--Step mode\n");
		jacobiMPIIterationStep_1D(n, &epsilon, step,
										b, m, w, &nTime1, &nTime2);
		printf("--Epsilon = %lf\n", epsilon);	
	}
	printf("--Result outputing...");
	char			*outDir = getOutDir(n, epsilon, b, step, outFile);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	//output result
	outMatrix1DtoF(m, n, outDir);
	//timer2 ends
	QueryPerformanceCounter(&nStopCounter);
	//get time
	nTime3 = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");

	printf("--(Time/s)Init=%lf, Computing=%lf, Data-saving=%lf, Total=%lf\n", 
		nTime1, nTime2, nTime3, nTime1 + nTime2 + nTime3);

	outLog(n, epsilon, step, b, nTime1, nTime2, nTime3, outFile, outDir);

	return;
}