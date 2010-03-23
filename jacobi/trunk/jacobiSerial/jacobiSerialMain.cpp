// jacobiSerial.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"

//********************************************************************************
int main(int argc, char* argv[])
{
	//paramenters
	int				n;
	double			epsilon;
	struct boundary	b;
	long			step;
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

	getchar();

	return 0;
}

