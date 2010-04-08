#include "stdafx.h"

//input parameters(float) from a .txt file named 'filename'
bool input_f(const char *filename, 
			 float *left, float *up, float *right, float *down,
			 int *x, int *y, 
			 float *epsilon, 
			 int *iteration,
			 char *outputfilename)
{
	bool t = false;

	FILE *fp;
	
	*epsilon = 0.0;
	*iteration = 0;

	if ((fp = fopen(filename, "r")) == NULL)
		printf("cannot open input file: %s ", filename);
	else
	{
		fscanf(fp, "%f %f %f %f\n", left, up, right, down);
		fscanf(fp, "%d %d\n", x, y);
		fscanf(fp, "%f\n", epsilon);
		fscanf(fp, "%d\n", iteration); 
		fscanf(fp, "%s\n", outputfilename);	
		t = true;
	}		

	fclose(fp);
	return t;
}


//input parameters(double) from a .txt file named 'filename'
bool input_d(const char *filename,
			 double *left, double *up, double *right, double *down,
			 int *x, int *y, 
			 double *epsilon, 
			 int *iteration, 
			 char *outputfilename)
{
	bool t = false;

	FILE *fp;
	
	*epsilon = 0.0;
	*iteration = 0;

	if ((fp = fopen(filename, "r")) == NULL)
		printf("cannot open input file: %s ", filename);
	else
	{
		fscanf(fp, "%lf %lf %lf %lf\n", left, up, right, down);
		fscanf(fp, "%d %d\n", x, y);
		fscanf(fp, "%lf\n", epsilon);
		fscanf(fp, "%d\n", iteration); 
		fscanf(fp, "%s\n", outputfilename);
		t = true;
	}		

	fclose(fp);
	return t;
}