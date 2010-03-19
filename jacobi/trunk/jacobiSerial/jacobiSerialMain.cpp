// jacobiSerial.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"

//********************************************************************************
//get input arguments. there is 3 ways -
//	way 1: read arguments from default .txt file named "input.txt";(no arguments)
//	way 2: read arguments from certain .txt file whose name is the first arguments 
//			after command;(1 argument)
//	way 3: read atguments directly from arguments.(9 arguments)
int input(int argc, char* argv[], 
		  int *n, double *epsilon, long int *step, struct boundary *b, 
		  char *outFile)
{
	if (argc == 1)
	{
		char		*filename = DEFAULT_INPUT_FILE;
		FILE		*fp;
		if ((fp = fopen(filename, "r")) == NULL)
		{
			printf("Cannot open defalut input file - %s.\n", filename);
			return 1;
		}
		else
		{
			fscanf(fp, "%lf %lf %lf %lf\n", 
				&(*b).left, &(*b).up, &(*b).right, &(*b).down);
			fscanf(fp, "%d\n", n);
			fscanf(fp, "%lf\n", epsilon);
			fscanf(fp, "%d\n", step);
			fscanf(fp, "%s\n", outFile);
			fclose(fp);
		}
	}
	else if (argc == 2)
	{
		FILE		*fp;
		if ((fp = fopen(argv[1], "r")) == NULL)
		{
			printf("Cannot open given input file - %s.\n", filename);
			return 2;
		}
		else
		{
			fscanf(fp, "%lf %lf %lf %lf\n", 
				&(*b).left, &(*b).up, &(*b).right, &(*b).down);
			fscanf(fp, "%d\n", n);
			fscanf(fp, "%lf\n", epsilon);
			fscanf(fp, "%d\n", step);
			fscanf(fp, "%s\n", outFile);
			fclose(fp);
		}
	}
	else if (argc == 9)
	{
		*n = atoi(argv[1]);
		*epsilon = atof(argv[2]);
		*step = stoi(argv[3]);
		(*b).left = atof(argv[4]);
		(*b).up = atof(argv[5]);
		(*b).right = atof(argv[6]);
		(*b).down = atof(argv[7]);
		strcpy(outFile, argv[8]);			
	}
	else
	{
		printf("The count of arguments is wrong. Check it. \n");
		return 3;
	}	

	return 0;
	//000000
}

//********************************************************************************
void jacobiSerial_1D(int n, double epsilon, 
					 long int step, struct boundary b, char *outFile)
{


}

//********************************************************************************
void jacobiSerialIterationEpsilon_2D(const int n, const double epsilon, 
									 int *step, const struct boundary b, 
									 double *u, double *w)
{
	double			*temp;
	//init data
	printf("Data initing...", n, epsilon);
	double			average = (b.left + b.up + b.right + b.down) / 4.0;
	for(int i = 0; i < n - 1; i++)
	{
		u[i * n + n] = b.left; 
		u[i * n + n - 1]= b.right; 
		u[i] = b.up; 
		u[(n - 1) * n + i + 1] = b.down;
		w[i * n + n] = b.left; 
		w[i * n + n - 1]= b.right; 
		w[i] = b.up; 
		w[(n - 1) * n + i + 1] = b.down;
	}
	for(int i = 1; i < n - 1; i++)
		for(int j = 1; j < n - 1; j++)
			u[i * n + j] = w[i * n + j] = average;
	printf("DONE.\n");
	//iteration
	printf("Computing...", n, epsilon);
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
				w[i * n + j] = (u[i * n + j - n] + u[i * n + j + n] 
							+ u[i * n + j - 1] + u[i * n + j + 1]) / 4.0;
				if (*step % 2 == 0)
				if(fabs(w[i * n + j] - u[i * n + j]) < epsilon) num++;
			}
		temp = u;
		u = w;
		w = temp;
	}
	printf("DONE.\n");
	return;

}

void jacobiSerialIterationStep_2D(const int n, const double *epsilon, 
								  int step, const struct boundary b, 
								  double *u, double *w)
{
}

//********************************************************************************
void jacobiSerial_2D(int n, double epsilon, 
					 long int step, struct boundary b, char *outFile)
{	
	//more paramenters
	double			*m = (double *)malloc(sizeof(double) * n * n);	
	double			*w = (double *)malloc(sizeof(double) * n * n);
	char			*resultFile;

	//timer
	LARGE_INTEGER	nFrequency, nStartCounter, nStopCounter;
	int				mul = 100000;
	double			nTime1, nTime2;
	QueryPerformanceFrequency(&nFrequency);

	//jacobi serial 2D solution
	if (epsilon == 0)
	{
		//timer starts
		QueryPerformanceCounter(&nStartCounter);

		jacobiSerialIterationEpsilon_2D(n, epsilon, &step, b, m, w);	

		//timer ends
		QueryPerformanceCounter(&nStopCounter);
	}
	else 
	{
		//timer starts
		QueryPerformanceCounter(&nStartCounter);

		jacobiSerialIterationStep_2D(n, &epsilon, b, step, m, w);	

		//timer ends
		QueryPerformanceCounter(&nStopCounter);
	}
	//get time
	nTime1 = (double)mul * (nStopCounter.QuadPart - nStartCounter.QuadPart) / nFrequency.QuadPart;
	nTime1 = nTime1 / mul;

	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	//output result
	outResultFilename = getOutResultFilename(n, epsilon, b, step, outFilename);
	ok = outputDoubleMatrixtoFile(matrix, n, outResultFilename);
	//timer2 ends
	QueryPerformanceCounter(&nStopCounter);
	//get time
	nTime2 = (double)mul * (nStopCounter.QuadPart - nStartCounter.QuadPart) / nFrequency.QuadPart;
	nTime2 = nTime1 / mul;

	return;
}

//********************************************************************************
int main(int argc, char* argv[])
{
	//paramenters
	int				n;
	double			epsilon;
	struct boundary	b;
	long int		step;
	char			*outFile = (char *)malloc(sizeof(char) * 128);

	//input
	int ok = input(argc, argv, &n, &epsilon, &step, &b, outFile);
	if (ok != 0)
	{
		getchar();
		return 0;
	}

	//jacobi serial 1D
	jacobiSerial_1D(n, epsilon, step, b, outFile);

	//jacobi serial 2D
	jacobiSerial_2D(n, epsilon, step, b, outFile);

	return 0;
}

