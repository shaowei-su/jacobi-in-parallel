#include "stdafx.h"


//*****************************************************************************
//Serial jacobi with FLOAT
//*****************************************************************************
//init with float
void serialInitMatrix_f(int x, int y, 
						float left, float up, float right, float down, 
						float *u, float *w)
{
	//init left & right boundary
	for(int i = 1; i < x - 1; i++)
	{
		u[i * y] = left;
		w[i * y] = left;
		u[i * y + y - 1] = right;
		w[i * y + y - 1] = right;
	}
	//init up & down boundary
	for(int i = 1; i < y - 1; i++)
	{
		u[i] = up;
		w[i] = up;
		u[x * y - y + i] = down;
		w[x * y - y + i] = down;
	}
	//init 4 corner.
	u[0] = up;
	u[y - 1] = right;
	u[x * y - 1] = down;
	u[x * y - y] = left;
	//init non-boundary
	float		average = ((x - 2) * (left + right) + (y - 2) * (up + down))
								/ ((x + y - 4) * 2);
	for(int i = 1; i < x - 1; i++)
		for(int j = 1; j < x - 1; j++)
		{
			u[i * y + j] = average;
			w[i * y + j] = average;
		}
	return;
}

//jacobi iteration process with float in epsilon mode
void serialJacobiIteration_epsilon_f(int x, int y, float epsilon, 
									 float *u, float *w, int *iteration)
{
	float		*temp;
	int			goal = (x - 2) * (y - 2);
	int			num = 0;
	int			step = 0;
	int			id;

	while(num < goal)
	{
		step++;
		num = 0;
		for(int i = 1;i < x - 1; i++)
			for(int j = 1;j < y - 1; j++)
			{
				id = i * y + j;
				w[id] = (float)((u[id - y] + u[id + y] 
								+ u[id - 1] + u[id + 1]) / 4.0);
				if(fabs(w[id] - u[id]) < epsilon) num++;
			}
		temp = u;
		u = w;
		w = temp;
	}
	*iteration = step;
	return;
}

//jacobi iteration process with float in iteration mode
void serialJacobiIteration_iteration_f(int x, int y, int iteration, 
									   float *u, float *w, float *epsilon)
{
	float		*temp;

	for(int k = 0; k < iteration; k++)
	{
		for(int i = 1;i < x - 1; i++)
			for(int j = 1;j < y - 1; j++)
				w[i * y + j] = (float)((u[i * y + j - y] + u[i * y + j + y] 
								+ u[i * y + j - 1] + u[i * y + j + 1]) / 4.0);
		temp = u;
		u = w;
		w = temp;
	}

	*epsilon = 0.0;
	for(int i = 1;i < x - 1; i++)
		for(int j = 1;j < y - 1; j++)
			if(fabs(w[i * y + j] - u[i * y + j]) > *epsilon) 
				*epsilon = fabs(w[i * y + j] - u[i * y + j]);

	return;			
}

//serial jacobi with float
void serialJacobi_f(float left, float up, float right, float down, 
					int x, int y, 
					float epsilon, 
					int iteration, 
					const char *outputfilename)
{
	//variables
	float		*u, *w;	
	u = (float *)malloc(sizeof(float) * x * y);
	w = (float *)malloc(sizeof(float) * x * y);

	//timer
	LARGE_INTEGER nFrequency;
	LARGE_INTEGER nStartCounter;
	LARGE_INTEGER nStopCounter;
	int			mul = 100000;
	double		nTime;
	QueryPerformanceFrequency(&nFrequency);

	//init
	QueryPerformanceCounter(&nStartCounter);//timer starts
	serialInitMatrix_f(x, y, left, up, right, down, u, w);
	QueryPerformanceCounter(&nStopCounter);//timer ends
	nTime = (double)mul * (nStopCounter.QuadPart - nStartCounter.QuadPart) 
													/ nFrequency.QuadPart;
	outputDouble(nTime / mul, " - Init time (second) = ", true);

	//jacobi iteration process in epsilon mode
	if (epsilon != 0 && iteration == 0)
	{
		outputString(" - Epsilon Mode ", true);
		outputFloat(epsilon, "\tEpsilon = ",true);
		outputString("\tComputing starts...", false);

		QueryPerformanceCounter(&nStartCounter);//timer starts
		serialJacobiIteration_epsilon_f(x, y, epsilon, u, w, &iteration);
		QueryPerformanceCounter(&nStopCounter);//timer ends
		nTime = (double)mul * (nStopCounter.QuadPart - nStartCounter.QuadPart) 
														/ nFrequency.QuadPart;
		outputDouble(nTime / mul, " - Epsilon Mode time (second) = ", true);

		outputString("Computing ends.", true);
		outputInt(iteration, "\tIteration count = ", true);
	}
	//jacobi iteration process in iteration mode
	else if (epsilon == 0 && iteration != 0)
	{
		outputString(" - Iteration Mode", true);
		outputInt(iteration, "\tIteration = ",true);
		outputString("\tComputing starts...", false);

		serialJacobiIteration_iteration_f(x, y, iteration, u, w, &epsilon);

		outputString("Computing ends.", true);
		outputFloat(epsilon, "\tEpsilon = ", true);
	}
	else
		outputString("Epsilon or Iteration is wrong.", true);

	//output matix
	//outputMatrix_f(u, x, y, "matrix from Serial :", false);
	QueryPerformanceCounter(&nStartCounter);//timer starts
	outputFloatMatrixtoFile(u, x, y, "matrix from Serial with Float :", outputfilename);
	QueryPerformanceCounter(&nStopCounter);//timer ends
	nTime = (double)mul * (nStopCounter.QuadPart - nStartCounter.QuadPart) 
													/ nFrequency.QuadPart;
	outputDouble(nTime / mul, " - Output to File time (second) = ", true);

	free(u);
	free(w);
}


//*****************************************************************************
//Serial jacobi with DOUBLE
//*****************************************************************************
//init with double
void serialInitMatrix_d(int x, int y, 
				  double left, double up, double right, double down, 
				  double *u, double *w)
{
	//init left & right boundary
	for(int i = 1; i < x - 1; i++)
	{
		u[i * y] = left;
		w[i * y] = left;
		u[i * y + y - 1] = right;
		w[i * y + y - 1] = right;
	}
	//init up & down boundary
	for(int i = 1; i < y - 1; i++)
	{
		u[i] = up;
		w[i] = up;
		u[x * y - y + i] = down;
		w[x * y - y + i] = down;
	}
	//init 4 corner.
	u[0] = up;
	u[y - 1] = right;
	u[x * y - 1] = down;
	u[x * y - y] = left;
	//init non-boundary
	double		average = ((x - 2) * (left + right) + (y - 2) * (up + down))
								/ ((x + y - 4) * 2);
	for(int i = 1; i < x - 1; i++)
		for(int j = 1; j < x - 1; j++)
		{
			u[i * y + j] = average;
			w[i * y + j] = average;
		}
	return;
}

//jacobi iteration process with double in epsilon mode
void serialJacobiIteration_epsilon_d(int x, int y, double epsilon, 
									 double *u, double *w, int *iteration)
{
	double		*temp;
	int			goal = (x - 2) * (y - 2);
	int			num = 0;
	int			step = 0;

	while(num < goal)
	{
		step++;
		num = 0;
		for(int i = 1;i < x - 1; i++)
			for(int j = 1;j < y - 1; j++)
			{
				w[i * y + j] = (u[i * y + j - y] + u[i * y + j + y] 
								+ u[i * y + j - 1] + u[i * y + j + 1]) / 4.0;
				if(fabs(w[i * y + j] - u[i * y + j]) < epsilon) num++;
			}
		temp = u; u = w; w = temp;
	}
	*iteration = step;
	return;
}

//jacobi iteration process with double in iteration mode
void serialJacobiIteration_iteration_d(int x, int y, int iteration, 
									   double *u, double *w, double *epsilon)
{
	double		*temp;

	for(int k = 0; k < iteration; k++)
	{
		for(int i = 1;i < x - 1; i++)
			for(int j = 1;j < y - 1; j++)
				w[i * y + j] = (u[i * y + j - y] + u[i * y + j + y] 
								+ u[i * y + j - 1] + u[i * y + j + 1]) / 4.0;
		temp = u; u = w; w = temp;
	}

	*epsilon = 0.0;
	for(int i = 1;i < x - 1; i++)
		for(int j = 1;j < y - 1; j++)
			if(fabs(w[i * y + j] - u[i * y + j]) > *epsilon) 
				*epsilon = fabs(w[i * y + j] - u[i * y + j]);

	return;			
}

//serial jacobi with double
void serialJacobi_d(double left, double up, double right, double down, 
					int x, int y, 
					double epsilon, 
					int iteration, 
					const char *outputfilename)
{
	//variables
	double		*u, *w;	
	u = (double *)malloc(sizeof(double) * x * y);
	w = (double *)malloc(sizeof(double) * x * y);

	//init
	serialInitMatrix_d(x, y, left, up, right, down, u, w);

	//jacobi iteration process in epsilon mode
	if (epsilon != 0 && iteration == 0)
	{
		outputString(" - Epsilon Mode", true);
		outputDouble(epsilon, "\tEpsilon = ",true);
		outputString("\tComputing starts...", false);

		serialJacobiIteration_epsilon_d(x, y, epsilon, u, w, &iteration);

		outputString("Computing ends.", true);
		outputInt(iteration, "\tIteration count = ", true);
	}
	//jacobi iteration process in iteration mode
	else if (epsilon == 0 && iteration != 0)
	{
		outputString(" - Iteration Mode", true);
		outputInt(iteration, "\tIteration = ",true);
		outputString("\tComputing starts...", false);

		serialJacobiIteration_iteration_d(x, y, iteration, u, w, &epsilon);

		outputString("Computing ends.", true);
		outputDouble(epsilon, "\tEpsilon = ", true);
	}
	//output matix
	//outputMatrix_f(u, x, y, "matrix from Serial :", false);
	outputDoubleMatrixtoFile(u, x, y, "matrix from Serial with Double :", outputfilename);

	free(u);
	free(w);
}