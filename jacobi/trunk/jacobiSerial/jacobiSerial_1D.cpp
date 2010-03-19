
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
									 long *step, const struct boundary b, 
									 double *m, double *w,
									 double *initTime, double *iterTime)
{
	double			*temp;
	
	//timer
	LARGE_INTEGER	nStartCounter, nStopCounter;

	//init data
	printf("--Data initing(%d, %lf).....", n, epsilon);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	matrix1DInit(n, b, m, w); 
	//timer ends
	QueryPerformanceCounter(&nStopCounter);
	*initTime = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");
	//iteration
	printf("--Computing(%d, %lf).....", n, epsilon);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	int				num = 0;
	int				goal = (n - 2) * (n - 2);
	*step = 0;
	while(num < goal)
	{
		(*step)++;
		num = 0;
		for(int i = 1;i < n - 1; i++)
			for(int j = 1; j < n - 1; j++)
			{
				w[i * n + j] = (m[i * n + j - n] + m[i * n + j + n] 
								+ m[i * n + j - 1] + m[i * n + j + 1]) / 4.0;
				if (*step % JUMP == 0)
					if(fabs(w[i * n + j] - m[i * n + j]) < epsilon) 
						num++;
			}
		temp = m; m = w; w = temp;
	}
	//timer ends
	QueryPerformanceCounter(&nStopCounter);
	*iterTime = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");

	return;
}


//********************************************************************************
double getEpsilon(const int n, double *m, double *w)
{
	double ep = 0;
	for(int i = 1; i < n - 1; i++)
		for(int j = 1; j < n - 1; j++)			
			if(fabs(w[i * n + j] - m[i * n + j]) < ep) 
				ep = fabs(w[i * n + j] - m[i * n + j]);
	return ep;
}

//********************************************************************************
void jacobiSerialIterationStep_1D(const int n, double *epsilon, 
								  const long step, const struct boundary b, 
								  double *m, double *w,
								  double *initTime, double *iterTime)
{
	double			*temp;
	
	//timer
	LARGE_INTEGER	nStartCounter, nStopCounter;

	//init data
	printf("--Data initing(%d, %d).....", n, step);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	matrix1DInit(n, b, m, w); 
	//timer ends
	QueryPerformanceCounter(&nStopCounter);
	*initTime = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");
	//iteration
	printf("--Computing(%d, %d).....", n, step);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	for( int k = 0; k < step; k++)
	{
		for(int i = 1;i < n - 1; i++)
			for(int j = 1; j < n - 1; j++)
				w[i * n + j] = (m[i * n + j - n] + m[i * n + j + n] 
								+ m[i * n + j - 1] + m[i * n + j + 1]) / 4.0;
		temp = m; m = w; w = temp;
	}
	*epsilon = getEpsilon(n, m, w);
	//timer ends
	QueryPerformanceCounter(&nStopCounter);
	*iterTime = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");

	return;
}

//********************************************************************************
void jacobiSerial_1D(int n, double epsilon, 
					 long step, struct boundary b, char *outFile)
{
	printf("Jacobi Serial 1D -\n");
	printf("--n=%d, e=%lf, step=%ld, LURD: %lf, %lf, %lf, %lf\n",
		n, epsilon, step, b.left, b.up, b.right, b.down);
	//more paramenters
	double			*m = (double *)malloc(sizeof(double) * n * n);
	double			*w = (double *)malloc(sizeof(double) * n * n);
	char			*outDir;

	//timer
	LARGE_INTEGER	nStartCounter, nStopCounter;
	double			nTime1, nTime2, nTime3;

	//jacobi serial 1D solution
	if (epsilon != 0)
	{
		printf("--Epsilon mode\n");
		jacobiSerialIterationEpsilon_1D(n, epsilon, &step, 
										b, m, w, &nTime1, &nTime2);
		printf("--Step = %ld", step);
	}
	else 
	{
		printf("--Step mode\n");
		jacobiSerialIterationStep_1D(n, &epsilon, step,
										b, m, w, &nTime1, &nTime2);
		printf("--Epsilon = %lf", epsilon);	
	}
	printf("--Result outputing...");
	outDir = getOutDir(n, epsilon, b, step, outFile);
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