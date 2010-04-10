// jacobiMPI.cpp : 定义控制台应用程序的入口点。
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

	b.averageValue = (b.left + b.up + b.right + b.down) / 4;
	
	//MPI_Init(&argc,&argv);

	//jacobi serial 1D
	jacobiMPI_1D(argc, argv, n, epsilon, step, b, outFile);

	//jacobi serial 2D
	//jacobiSerial_2D(n, epsilon, step, b, outFile);

	//MPI_Finalize();

	getchar();

	return 0;
}


