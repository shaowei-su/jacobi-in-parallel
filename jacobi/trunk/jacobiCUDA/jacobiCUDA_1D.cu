#include "stdafx.h"

#include <cuda_runtime.h>
#include <cutil_inline.h>


//*****************************************************************************
//init matrix non-boundarys
__global__ void kernelInitMatrixInner(const int n, const double value,
									  double *m, double *w)
{
	int id = blockIdx.x * blockDim.x +threadIdx.x;
	if (id < n * n)
	{
		m[id] = value;
		w[id] = value;
	}
}

//*****************************************************************************
//init matrix boundarys
__global__ void kernelInitMatrixBoundary(const int n, const boundary b,
										 double *m, double *w)
{
	int id = blockIdx.x * blockDim.x +threadIdx.x;
	if (id < n - 1)
	{
		int leftId = id * n + n;
		m[leftId] = b.left; w[leftId] = b.left;
		int upId = id;
		m[upId] = b.up;	w[upId] = b.up;		
		int rightId = leftId - 1;
		m[rightId] = b.right; w[rightId] = b.right;
		int downId = n * n - n + 1 + upId;
		m[downId] = b.down;	w[downId] = b.down;
	}
}

//*****************************************************************************
//init matrix
void initMatrix(const int n, const struct boundary b, double *d_m, double *d_w)
{
	int				blockSize = 256;
	dim3			block(blockSize);
	dim3			matrixBoundaryGrid(getQuotient(n - 1, blockSize));
	dim3			matrixInnerGrid(getQuotient(n * n, blockSize));

	kernelInitMatrixInner
		<<<matrixInnerGrid, block>>>(n, b.averageValue, d_m, d_w);
	kernelInitMatrixBoundary
		<<<matrixBoundaryGrid, block>>>(n, b, d_m, d_w);
}

//*****************************************************************************
//kernel of jacobi iteration process with double in iteration mode
__global__ void kernelJacobiIteration(const int n, double *m, double *w)
{
	int id = blockIdx.x * blockDim.x +threadIdx.x;
	int row = id / (n - 1);
	int column = id - row * (n - 1);
	int location = (row + 1) * (column + 1);
	w[location] = (m[location - 1] + m[location - n] + m[location + 1] 
										+ m[location + n]) / 4.0;
}

//*****************************************************************************
//kernel of getting epsilon between d_m and d_w
__global__ void kernelGetEpsilon(const int n, 
								 double *m, double *w, 
								 double *ep)
{
	int id = blockIdx.x * blockDim.x +threadIdx.x;
	int row = id / (n - 1);
	int column = id - row * (n - 1);
	int location = (row + 1) * (column + 1);
	ep[id] = fabs(m[location] - w[location]);	
}

//kernel of max of a num group
__global__ void kernelGetMax(const int count, double *ep)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	if (id < count)
		if (ep[id] < ep[id + count]) ep[id] = ep[id + count];
}

//*****************************************************************************
//get epsilon 
double getEpsilon(const int n, double *d_m, double *d_w)
{
	double			epsilon;
	double			*d_ep;
    //cutilSafeCall
	(cudaMalloc((void**) &d_ep, (n - 1) * (n - 1) * sizeof(double)));
	// setup execution parameters	
	int				blockSize = 256;
    dim3			block(blockSize);
    dim3			matrixInnerGrid(getQuotient((n - 1) * (n - 1), blockSize));
	dim3			getMaxGrid(getQuotient((n - 1) * (n - 1) / 2, blockSize));
	kernelGetEpsilon<<<matrixInnerGrid, block>>>(n, d_m, d_w, d_ep);
	//check if kernel execution generated and error
	//cutilCheckMsg("Kernel execution failed");
	//get max of d_ep
	int			count = (int)((n - 1) * (n - 1) / 2);
	while(count > 1)
	{
		kernelGetMax<<<getMaxGrid, block>>>(count, d_ep);
		//check if kernel execution generated and error
		//cutilCheckMsg("Kernel execution failed");
		count = (count + 1) / 2;
	}

	double			lastD_ep;
	//cutilSafeCall
		(cudaMemcpy(&lastD_ep, d_ep + (n - 1) * (n - 1) - 1, 
									sizeof(double), cudaMemcpyDeviceToHost));
	//cutilSafeCall
		(cudaMemcpy(&epsilon, d_ep, sizeof(double), cudaMemcpyDeviceToHost));
	//cutilSafeCall
		(cudaFree(d_ep));

	if (epsilon < lastD_ep) epsilon = lastD_ep;

	return epsilon;
}

//*****************************************************************************
void jacobiCUDAIterationEpsilon_1D(const int n, const double epsilon, 
								   long *step, const struct boundary b, 
								   double *d_m, double *d_w,
								   double *initTime, double *iterTime)
{
	double			*temp;
	//timer
	LARGE_INTEGER	nStartCounter, nStopCounter;	
	//init data
	printf("--Data initing(%d, %lf).....", n, epsilon);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	//init
	initMatrix(n, b, d_m, d_w);

//*****************************************************************************
	double			*m = (double *)malloc(sizeof(double) * n * n);
	//cutilSafeCall
	(cudaMemcpy(m, d_m, matrixMemSize, cudaMemcpyDeviceToHost));
	printf("\n");
	for(int i = 0; i < 10; i++)
		for(int j = 0; j < 10; j++)
			printf("%.2lf  ", m[i * n + j]);
	
//*****************************************************************************

	//timer ends
	QueryPerformanceCounter(&nStopCounter);
	*initTime = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");
	//iteration
	printf("--Computing(%d, %lf).....", n, epsilon);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	*step = 0;
	double			epsilonTemp = epsilon + 1;
	// setup execution parameters	
	int				blockSize = 256;
    dim3			block(blockSize);
    dim3			matrixInnerGrid(getQuotient((n - 1) * (n - 1), blockSize));
	//iteration
    while (epsilonTemp > epsilon)
	{
		(*step)++;
		printf("%d\n", *step);
		//execute the kernel
		kernelJacobiIteration<<<matrixInnerGrid, block>>>(n, d_m, d_w);
		//check if kernel execution generated and error
		//cutilCheckMsg("Kernel execution failed");
		//epsilon
		if (*step % JUMP == 0)
			epsilonTemp = getEpsilon(n, d_m, d_w);
		temp = d_m; d_m = d_w; d_w = temp;
    }
	//timer ends
	QueryPerformanceCounter(&nStopCounter);
	*iterTime = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");
	return;
}

//********************************************************************************
void jacobiCUDAIterationStep_1D(const int n, double *epsilon, 
								const long step, const struct boundary b,
								double *d_m, double *d_w,
								double *initTime, double *iterTime)
{
	double			*temp;
	//timer
	LARGE_INTEGER	nStartCounter, nStopCounter;	
	//init data
	printf("--Data initing(%d, %lf).....", n, epsilon);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	//init
	initMatrix(n, b, d_m, d_w);
	//timer ends
	QueryPerformanceCounter(&nStopCounter);
	*initTime = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");
	//iteration
	printf("--Computing(%d, %lf).....", n, epsilon);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	// setup execution parameters	
	int				blockSize = 256;
    dim3			block(blockSize);
    dim3			matrixInnerGrid(getQuotient((n - 1) * (n - 1), blockSize));
	for(int i = 0; i < step; i++)
	{
		//execute the kernel
		kernelJacobiIteration<<<matrixInnerGrid, block>>>(n, d_m, d_w);
		temp = d_m; d_m = d_w; d_w = temp;		
	}
	*epsilon = getEpsilon(n, d_m, d_w);
	//timer ends
	QueryPerformanceCounter(&nStopCounter);
	*iterTime = getCostTime(nStartCounter, nStopCounter);
	printf("Done.\n");
	return;
}

//********************************************************************************
void jacobiCUDA_1D(int argc, char** argv,
				   int n, double epsilon, 
				   long step, struct boundary b, char *outFile)
{
	printf("Jacobi CUDA 1D -\n");
	printf("--n=%d, e=%lf, step=%ld\n--LURD: %lf, %lf, %lf, %lf\n",
		n, epsilon, step, b.left, b.up, b.right, b.down);

	//init cuda device
    if (cutCheckCmdLineFlag(argc, (const char**)argv, "device"))
        cutilDeviceInit(argc, argv);
    else
        cudaSetDevice(cutGetMaxGflopsDeviceId());
	
    //allocate devce memory for matrices u
    unsigned int	matrixSize = n * n;
    unsigned int	matrixMemSize = sizeof(double) * matrixSize;

    double			*d_m;
	cudaMalloc((void**) &d_m, matrixMemSize);
    //cutilSafeCall(cudaMalloc((void**) &d_m, matrixMemSize));
    double			*d_w;
    //cutilSafeCall
		(cudaMalloc((void**) &d_w, matrixMemSize));

	//timer
	LARGE_INTEGER	nStartCounter, nStopCounter;
	double			nTime1, nTime2, nTime3;

	//jacobi serial 1D solution
	if (epsilon != 0)
	{
		printf("--Epsilon mode\n");
		jacobiCUDAIterationEpsilon_1D(n, epsilon, &step, 
										b, d_m, d_w, &nTime1, &nTime2);
		printf("--Step = %ld\n", step);
	}
	else 
	{
		printf("--Step mode\n");
		jacobiCUDAIterationStep_1D(n, &epsilon, step,
										b, d_m, d_w, &nTime1, &nTime2);
		printf("--Epsilon = %lf\n", epsilon);	
	}
	printf("--Result outputing...");

	char			*outDir = getOutDir(n, epsilon, b, step, outFile);
	//timer starts
	QueryPerformanceCounter(&nStartCounter);
	double			*m = (double *)malloc(sizeof(double) * n * n);
	//cutilSafeCall
		(cudaMemcpy(m, d_m, matrixMemSize, cudaMemcpyDeviceToHost));
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