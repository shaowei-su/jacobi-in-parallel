
#include "stdafx.h"

//********************************************************************************
void matrix1DInit(const int n, const struct boundary b, double *m, double *w)
{	
	double average = (b.left + b.up + b.right + b.down) / 4.0;
	for(int i = 0; i < n - 1; i++)
	{
		m[i * n + n]			= w[i * n + n]			= b.left; 
		m[i * n + n - 1]		= w[i * n + n - 1]		= b.right; 
		m[i]					= w[i]					= b.up; 
		m[n * n - n + i + 1]	= w[n * n - n + i + 1]	= b.down;
	}
	for(int i = 1; i < n - 1; i++)
		for(int j = 1; j < n - 1; j++)
			m[i * n + j] = w[i * n + j] = average;
	return;
}

//********************************************************************************
void jacobiSerialIterationEpsilon_1D(const int n, const double epsilon, 
									 int *step, const struct boundary b, 
									 double *m, double *w,
									 double *initTime, double *iterTime)
{
	double			*temp;
	
	//timer
	LARGE_INTEGER	nStartCounter, nStopCounter;
	//init data
	matrix1DInit(n, b, m, w); 
	printf("--Data initing(%d, %lf).....", n, epsilon);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	//timer ends
	QueryPerformanceCounter(&nStopCounter);
	*initTime = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");
	//iteration
	printf("--Computing(%d, %lf).....", n, epsilon);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	int				num;
	int				goal = (n - 2) * (n - 2);
	*step = 0;
	while(num < goal)
	{
		(*step)++;
		num = 0;
		for(int i = 1;i < n - 1; i++)
			for(int j = 1; j < n - 1; j++)
			{
				w[i * n + j] = (u[i * n + j - n] + u[i * n + j + n] 
								+ u[i * n + j - 1] + u[i * n + j + 1]) / 4.0;
				//if (*step % 2 == 0)
				if(fabs(w[i * n + j] - u[i * n + j]) < epsilon) num++;
			}
		temp = u; u = w; w = temp;
	}
	//timer ends
	QueryPerformanceCounter(&nStopCounter);
	*iterTime = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");

	return;
}

//********************************************************************************
void jacobiSerialIterationStep_1D(const int n, double *epsilon, 
								  const int step, const struct boundary b, 
								  double *u, double *w)
{
	double			*temp;
	
	//timer
	LARGE_INTEGER	nStartCounter, nStopCounter;

	//init data
	printf("--Data initing(%d, %lf).....", n, epsilon);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	double			average = (b.left + b.up + b.right + b.down) / 4.0;
	for(int i = 0; i < n - 1; i++)
	{
		u[i * n + n]			= w[i * n + n]			= b.left; 
		u[i * n + n - 1]		= w[i * n + n - 1]		= b.right; 
		u[i]					= w[i]					= b.up; 
		u[n * n - n + i + 1]	= w[n * n - n + i + 1]	= b.down;
	}

	for(int i = 1; i < n - 1; i++)
		for(int j = 1; j < n - 1; j++)
			u[i * n + j] = w[i * n + j] = average;
	//timer ends
	QueryPerformanceCounter(&nStopCounter);
	*initTime = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");
	//iteration
	printf("--Computing(%d, %lf).....", n, epsilon);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	for( int i = 0; i < step; i++)
	{
		for(int i = 1;i < n - 1; i++)
			for(int j = 1; j < n - 1; j++)
				w[i * n + j] = (u[i * n + j - n] + u[i * n + j + n] 
								+ u[i * n + j - 1] + u[i * n + j + 1]) / 4.0;
		temp = u; u = w; w = temp;
	}
	//timer ends
	QueryPerformanceCounter(&nStopCounter);
	*iterTime = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");

	return;
}

//********************************************************************************
void jacobiSerial_1D(int n, double epsilon, 
					 long int step, struct boundary b, char *outFile)
{
	//more paramenters
	double			*m = (double *)malloc(sizeof(double) * n * n);
	double			*w = (double *)malloc(sizeof(double) * n * n);
	char			*resultFile;

	//timer
	LARGE_INTEGER	nStartCounter, nStopCounter;
	double			nTime1, nTime2, nTime3;

	//jacobi serial 1D solution
	if (epsilon == 0)
	{
		printf("Jacobi Serial 1D - Epsilon mode\n");
		jacobiSerialIterationEpsilon_2D(n, epsilon, &step, 
										b, m, w, &nTime1, &nTime2);	
	}
	else 
	{
		printf("Jacobi Serial 1D - Step mode\n");
		jacobiSerialIterationStep_2D(n, &epsilon, step,
										b, m, w, &nTime1, &nTime2);	
	}

	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	//output result
	outResultFilename = getOutResultFilename(n, epsilon, b, step, outFilename);
	ok = outputDoubleMatrixtoFile(matrix, n, outResultFilename);
	//timer2 ends
	QueryPerformanceCounter(&nStopCounter);
	//get time
	nTime3 = getCostTime(nStartCounter, nStopCounter);

	return;
}