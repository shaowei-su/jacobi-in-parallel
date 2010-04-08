#include "stdafx.h"

int main(int argc, char** argv)
{
	//float
	float		left_f, up_f, right_f, down_f;
	int			x, y;
	float		epsilon_f;
	int			iteration;
	char		*outputfilename = (char*)malloc(sizeof(char)*32);
	char		*outputfilename2 = (char*)malloc(sizeof(char)*32);


	//double
	double		left_d, up_d, right_d, down_d;
	double		epsilon_d;

	//timer
	LARGE_INTEGER nFrequency;
	LARGE_INTEGER nStartCounter;
	LARGE_INTEGER nStopCounter;
	int			mul = 100000;
	double		nTime;
	QueryPerformanceFrequency(&nFrequency);

	/////////////////////////////////////////
	////FLOAT
	/////////////////////////////////////////
	////input with FLOAT
	//QueryPerformanceCounter(&nStartCounter);//timer starts
	//if (input_f("input_f.txt", &left_f, &up_f, &right_f, &down_f, 
	//				&x, &y, &epsilon_f, &iteration, outputfilename))
	//	test_output_f("input_f.txt", left_f, up_f, right_f, down_f, 
	//					x, y, epsilon_f, iteration, outputfilename);
	//else
	//	outputString("Input error @ input_f.txt.", true);	
	//QueryPerformanceCounter(&nStopCounter);//timer ends
	//nTime = (double)mul * (nStopCounter.QuadPart - nStartCounter.QuadPart) 
	//												/ nFrequency.QuadPart;
	//outputDouble(nTime / mul, " - Input FLOAT time (second) = ", true);
	//outputString("\n", true);

	////serial jacobi with float
	//strcpy(outputfilename2, outputfilename);
	//strcat(outputfilename2, "_Serial.txt");
	//QueryPerformanceCounter(&nStartCounter);//timer starts
	//outputString("Serial Jacobi with FLOAT", true);
	//serialJacobi_f(left_f, up_f, right_f, down_f, 
	//				x, y, epsilon_f, iteration, outputfilename2);
	//QueryPerformanceCounter(&nStopCounter);//timer ends
	//nTime = (double)mul * (nStopCounter.QuadPart - nStartCounter.QuadPart) 
	//												/ nFrequency.QuadPart;
	//outputDouble(nTime / mul, " - Serial FLOAT time (second) = ", true);
	//outputString("\n", true);

	////cuda jacobi with float
	//strcpy(outputfilename2, outputfilename);
	//strcat(outputfilename2, "_cuda.txt");
	//QueryPerformanceCounter(&nStartCounter);//timer starts
	//outputString("Cuda Jacobi with FLOAT", true);
	//cudaJacobi_f(argc, argv, left_f, up_f, right_f, down_f, 
	//					x, y, epsilon_f, iteration, outputfilename2);
	//QueryPerformanceCounter(&nStopCounter);//timer ends
	//nTime = (double)mul * (nStopCounter.QuadPart - nStartCounter.QuadPart) 
	//												/ nFrequency.QuadPart;
	//outputDouble(nTime / mul, " - Cuda FLOAT time (second) = ", true);
	//outputString("\n", true);

	///////////////////////////////////////
	//DOUBLE
	///////////////////////////////////////
	//input with DOUBLE
	QueryPerformanceCounter(&nStartCounter);//timer starts
	if (input_d("input_d.txt", &left_d, &up_d, &right_d, &down_d, 
					&x, &y, &epsilon_d, &iteration, outputfilename))
		test_output_d("input_d.txt", left_d, up_d, right_d, down_d,	
						x, y, epsilon_d, iteration, outputfilename);
	else
		outputString("Input error @ input_d.txt.", true);
	QueryPerformanceCounter(&nStopCounter);//timer ends
	nTime = (double)mul * (nStopCounter.QuadPart - nStartCounter.QuadPart) 
													/ nFrequency.QuadPart;
	outputDouble(nTime / mul, " - Input DOUBLE time (second) = ", true);
	outputString("\n", true);

	////serial jacobi with double
	//strcpy(outputfilename2, outputfilename);
	//strcat(outputfilename2, "_Serial.txt");
	//QueryPerformanceCounter(&nStartCounter);//timer starts
	//outputString("Serial Jacobi with DOUBLE", true);
	//serialJacobi_d(left_d, up_d, right_d, down_d, 
	//				x, y, epsilon_d, iteration, outputfilename2);
	//QueryPerformanceCounter(&nStopCounter);//timer ends
	//nTime = (double)mul * (nStopCounter.QuadPart - nStartCounter.QuadPart) 
	//												/ nFrequency.QuadPart;
	//outputDouble(nTime / mul, " - Serial DOUBLE time (second) = ", true);
	//outputString("\n", true);

	//cuda jacobi with double
	strcpy(outputfilename2, outputfilename);
	strcat(outputfilename2, "_cuda.txt");
	QueryPerformanceCounter(&nStartCounter);//timer starts
	outputString("Cuda Jacobi with DOUBLE", true);
	cudaJacobi_d(argc, argv, left_d, up_d, right_d, down_d, 
						x, y, epsilon_d, iteration, outputfilename2);
	QueryPerformanceCounter(&nStopCounter);//timer ends
	nTime = (double)mul * (nStopCounter.QuadPart - nStartCounter.QuadPart) 
													/ nFrequency.QuadPart;
	outputDouble(nTime / mul, " - Cuda DOUBLE time (second) = ", true);
	outputString("\n", true);
	outputTimetoFile(nTime / mul, "DOUBLE time (secong) = ", outputfilename2);	

	getchar();

	return 0;
}