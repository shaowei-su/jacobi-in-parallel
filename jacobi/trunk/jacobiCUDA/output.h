#include "stdafx.h"


void outputTimetoFile(double time, char* name, const char* filename);

//output a float matrix. the size of matrix is x * y
void outputFloatMatrix(float* m, int x, int y, char* name, bool flag);

//output a matrix with float into a file. the size of matrix is x * y
void outputFloatMatrixtoFile(float* m, int x, int y, char* name, const char* filename);

//output a double matrix. the size of matrix is x * y
void outputDoubleMatrix(double* m, int x, int y, char* name, bool flag);

//output a matrix with double into a file. the size of matrix is x * y
void outputDoubleMatrixtoFile(double* m, int x, int y, char* name, const char* filename);

//output a int.
void outputInt(int value, char *name, bool flag);

//output a int into a file.
void outputInttoFile(int value, char *name, bool flag, const char* filename);

//output a float.
void outputFloat(float value, char *name, bool flag);

//output a float into a file.
void outputFloattoFile(int value, char *name, bool flag, const char* filename);

//output a double.
void outputDouble(double value, char *name, bool flag);

//output a double into a file.
void outputDoubletoFile(int value, char *name, bool flag, const char* filename);

//output a string.
void outputString(char *name, bool flag);

//output a double into a file.
void outputStringtoFile(char *name, bool flag, const char* filename);

//test output
void test_output_f(const char *filename,
				   float left, float up, float right, float down,
				   int x, int y, 
				   float epsilon, 
				   int iteration, 
				   char *outputfilename);


//test output
void test_output_d(const char *filename,
				   double left, double up, double right, double down,
				   int x, int y, 
				   double epsilon, 
				   int iteration, 
				   char *outputfilename);