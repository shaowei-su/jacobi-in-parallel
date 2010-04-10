
#include "stdafx.h"

//********************************************************************************
void matrix2DInit(const int n, const struct boundary b, double **m, double **w)
{
	#pragma omp parallel for
	for(int i = 0; i < n - 1; i++)
	{
		m[i + 1][0]			= w[i + 1][0]		= b.left; 
		m[i][n - 1]			= w[i][n - 1]		= b.right; 
		m[0][i]				= w[0][i]			= b.up; 
		m[n - 1][i + 1]		= w[n - 1][i + 1]	= b.down;
	}
	for(int i = 1; i < n - 1; i++)
		#pragma omp parallel for
		for(int j = 1; j < n - 1; j++)
			m[i][j] = w[i][j] = b.averageValue;
	return;
}

//********************************************************************************
void jacobiSerialIterationEpsilon_2D(const int n, const double epsilon, 
									 long *step, const struct boundary b, 
									 double **m, double **w,
									 double *initTime, double *iterTime)
{
	double			**temp;
	
	//timer
	LARGE_INTEGER	nStartCounter, nStopCounter;

	//init data
	printf("--Data initing(%d, %lf).....", n, epsilon);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	matrix2DInit(n, b, m, w); 
	//timer ends
	QueryPerformanceCounter(&nStopCounter);
	*initTime = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");
	//iteration
	printf("--Computing(%d, %lf).....", n, epsilon);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	int				count;
	int				countGoal = (n - 2) * (n - 2);
	*step = 0;
	while(count < countGoal)
	{
		(*step)++;
		count = 0;
		for(int i = 1;i < n - 1; i++)
			#pragma omp parallel for reduction (+:count)
			for(int j = 1; j < n - 1; j++)
			{
				w[i][j] = (m[i - 1][j] + m[i + 1][j] 
								+ m[i][j - 1] + m[i][j + 1]) / 4.0;
				if (*step % JUMP == 0)
					if(fabs(w[i][j] - m[i][j]) < epsilon) 
						count++;
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
double getEpsilon_2D(const int n, double **m, double **w)
{
	double ep = 0;
	for(int i = 1; i < n - 1; i++)
		#pragma omp parallel for
		for(int j = 1; j < n - 1; j++)
			//#pragma omp critical			
			if(fabs(w[i][j] - m[i][j]) > ep) 
				ep = fabs(w[i][j] - m[i][j]);
	return ep;
}

//********************************************************************************
void jacobiSerialIterationStep_2D(const int n, double *epsilon, 
								  const long step, const struct boundary b, 
								  double **m, double **w,
								  double *initTime, double *iterTime)
{
	double			**temp;
	
	//timer
	LARGE_INTEGER	nStartCounter, nStopCounter;

	//init data
	printf("--Data initing(%d, %d).....", n, step);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	matrix2DInit(n, b, m, w); 
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
			#pragma omp parallel for
			for(int j = 1; j < n - 1; j++)
				w[i][j] = (m[i - 1][j] + m[i + 1][j] 
								+ m[i][j - 1] + m[i][j + 1]) / 4.0;
		temp = m; m = w; w = temp;
	}
	*epsilon = getEpsilon_2D(n, m, w);
	//timer ends
	QueryPerformanceCounter(&nStopCounter);
	*iterTime = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");

	return;
}

//********************************************************************************
void jacobiSerial_2D(int n, double epsilon, 
					 long step, struct boundary b, char *outFile)
{
	printf("Jacobi Serial 2D -\n");
	printf("--n=%d, e=%lf, step=%ld\nLURD: %lf, %lf, %lf, %lf\n",
		n, epsilon, step, b.left, b.up, b.right, b.down);
	//more paramenters
	double			**m = (double **)malloc(sizeof(double *) * n);
	double			**w = (double **)malloc(sizeof(double *) * n);
	for(int i = 0; i < n; i++)
	{
		m[i] = (double *)malloc(sizeof(double) * n);
		w[i] = (double *)malloc(sizeof(double) * n);
	}

	//timer
	LARGE_INTEGER	nStartCounter, nStopCounter;
	double			nTime1, nTime2, nTime3;

	//jacobi serial 1D solution
	if (epsilon != 0)
	{
		printf("--Epsilon mode\n");
		jacobiSerialIterationEpsilon_2D(n, epsilon, &step, 
										b, m, w, &nTime1, &nTime2);
		printf("--Step = %ld", step);
	}
	else 
	{
		printf("--Step mode\n");
		jacobiSerialIterationStep_2D(n, &epsilon, step,
										b, m, w, &nTime1, &nTime2);
		printf("--Epsilon = %lf", epsilon);	
	}
	printf("--Result outputing...");
	char			*outDir = getOutDir(n, epsilon, b, step, outFile);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	//output result
	outMatrix2DtoF(m, n, outDir);
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