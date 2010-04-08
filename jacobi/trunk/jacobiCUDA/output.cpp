#include "stdafx.h"

void outputTimetoFile(double time, char* name, const char* filename)
{
	FILE *fp;
	if((fp = fopen(filename, "at")) == NULL)
        printf("cannot open %s \n", filename);
	else
	{
		fprintf(fp, "\n***** %s %lf\n", name, time);
	}
	fclose(fp);
	return;
}

//output a matrix with float. the size of matrix is x * y
void outputFloatMatrix(float* m, int x, int y, char* name, bool flag)
{
	printf("Value of Matrix %s :\n", name);
	for(int i = 0; i < 10; i++)		
	{
		for(int j = 0; j < 10; j++)			
			printf("%.2f ", m[i * y + j]);
		printf("\n");
		if (flag) getchar();
	}
	printf("\n");
}

//output a matrix with float into a file. the size of matrix is x * y
void outputFloatMatrixtoFile(float* m, int x, int y, char* name, const char* filename)
{
	FILE *fp;
	if((fp = fopen(filename, "at")) == NULL)
        printf("cannot open %s ", filename);
	else
	{
		fprintf(fp, "Value of Matrix %s \n", name);
		for(int i = 0; i < x; i++)		
		{
			for(int j = 0; j < y; j++)			
				fprintf(fp, "%f ", m[i * y + j]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return;
}

//output a matrix with DOUBLE. the size of matrix is x * y
void outputDoubleMatrix(double* m, int x, int y, char* name, bool flag)
{
	printf("Value of Matrix %s :\n", name);
	for(int i = 0; i < 10; i++)		
	{
		for(int j = 0; j < 10; j++)			
			printf("%5lf ", m[i * y + j]);
		printf("\n");
		if (flag) getchar();
	}
	printf("\n");
}

//output a matrix with double into a file. the size of matrix is x * y
void outputDoubleMatrixtoFile(double* m, int x, int y, char* name, const char* filename)
{
	FILE *fp;
	if((fp = fopen(filename, "at")) == NULL)
        printf("cannot open %s ", filename);
	else
	{
		fprintf(fp, "Value of Matrix %s \n", name);
		for(int i = 0; i < x; i++)		
		{
			//fprintf(fp, "%d -  ", i);
			for(int j = 0; j < y; j++)			
				fprintf(fp, "%lf ", m[i * y + j]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return;
}

//output a int.
void outputInt(int value, char *name, bool flag)
{
	printf("%s %d ", name, value);
	if (flag) printf("\n");
	return;
}

//output a int into a file.
void outputInttoFile(int value, char *name, bool flag, const char* filename)
{
	FILE *fp;
	if((fp = fopen(filename, "at")) == NULL)
        printf("cannot open %s ", filename);
	else
	{
		fprintf(fp, "%s %d ", name, value);
		if (flag) fprintf(fp, "\n");
	}
	fclose(fp);
	return;
}

//output a float.
void outputFloat(float value, char *name, bool flag)
{
	printf("%s %f ", name, value);
	if (flag) printf("\n");
	return;
}

//output a float into a file.
void outputFloattoFile(int value, char *name, bool flag, const char* filename)
{
	FILE *fp;
	if((fp = fopen(filename, "at")) == NULL)
        printf("cannot open %s ", filename);
	else
	{
		fprintf(fp, "%s %f ", name, value);
		if (flag) fprintf(fp, "\n");
	}
	fclose(fp);
	return;
}

//output a double.
void outputDouble(double value, char *name, bool flag)
{
	printf("%s %lf ", name, value);
	if (flag) printf("\n");
	return;
}

//output a double into a file.
void outputDoubletoFile(int value, char *name, bool flag, const char* filename)
{
	FILE *fp;
	if((fp = fopen(filename, "at")) == NULL)
        printf("cannot open %s ", filename);
	else
	{
		fprintf(fp, "%s %lf ", name, value);
		if (flag) fprintf(fp, "\n");
	}
	fclose(fp);
	return;
}

//output a string.
void outputString(char *name, bool flag)
{
	printf("%s", name);
	if (flag) printf("\n");
	return;
}

//output a double into a file.
void outputStringtoFile(char *name, bool flag, const char* filename)
{
	FILE *fp;
	if((fp = fopen(filename, "at")) == NULL)
        printf("cannot open %s ", filename);
	else
	{
		fprintf(fp, "%s", name);
		if (flag) fprintf(fp, "\n");
	}
	fclose(fp);
	return;
}

//test output
void test_output_f(const char *filename, 
				   float left, float up, float right, float down,
				   int x, int y, 
				   float epsilon, 
				   int iteration, 
				   char *outputfilename)
{
	printf("from file: %s \n", filename);
	printf("\tleft = %f  \n\tup = %f  \n\tright = %f  \n\tdown = %f \n", 
		left, up, right, down);
	printf("\tx = %d  y = %d \n", x, y);
	printf("\tepsilon = %f  \n", epsilon);
	printf("\titeration = %d  \n", iteration);
	printf("\toutput filename: %s \n", outputfilename);
	printf("\n");
}


//test output
void test_output_d(const char *filename, 
				   double left, double up, double right, double down,
				   int x, int y, 
				   double epsilon, 
				   int iteration, 
				   char *outputfilename)
{
	printf("from file: %s \n", filename);
	printf("\tleft = %lf  \n\tup = %lf  \n\tright = %lf  \n\tdown = %lf \n", 
		left, up, right, down);
	printf("\tx = %d  y = %d \n", x, y);
	printf("\tepsilon = %7.15lf  \n", epsilon);
	printf("\titeration = %d  \n", iteration);
	printf("\toutput filename: %s \n", outputfilename);
	printf("\n");
}