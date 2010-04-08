#include "stdafx.h"

#include <cuda_runtime.h>
#include <cutil_inline.h>

int carry(int dividend, int divsor)
{
	int quotient = dividend / divsor;
	return quotient * divsor == dividend ? quotient : quotient + 1;
}

//*****************************************************************************
//Cuda jacobi with FLOAT
//*****************************************************************************
//init non-boundarys
__global__ void 
kernelInitInner_f(float* u, float* w, float average, int x, int y)
{
	int id = blockIdx.x * blockDim.x +threadIdx.x;
	if (id < x * y)
	{
		u[id] = average; w[id] = average;
	}
}

//init up & down boundarys
__global__ void kernelInitUpDownBoundary_f(float* u, float* w, 
									       float up, float down, int x, int y)
{
	int id = blockIdx.x * blockDim.x +threadIdx.x;
		if (id < y - 1)
	{
		int upLoca = id;
		u[upLoca] = up;	w[upLoca] = up;
		int downLoca = x * y - y + upLoca;
		u[downLoca] = down;	w[downLoca] = down;
	}
}

//init left & right boundarys
__global__ void kernelInitLeftRightBoundary_f(float* u, float* w, 
										float left, float right, int x, int y)
{
	int id = blockIdx.x * blockDim.x +threadIdx.x;
	if (id < x - 1)
	{
		int leftLoca = y + id * y;
		u[leftLoca] = left;	w[leftLoca] = left;
		int rightLoca = leftLoca - 1;
		u[rightLoca] = right; w[rightLoca] = right;
	}
}

//init with float
void cudaInitMatrix_f(int x, int y, 
				      float left, float up, float right, float down, 
				      float *d_u, float *d_w)
{
	int			size_block = 256;
	dim3		block(size_block);
	dim3		gridLine(carry(y - 1, size_block));
	dim3		gridColumn(carry(x - 1, size_block));
	dim3		gridInner(carry(x * y, size_block));

	float		average = ((x - 2) * (left + right) + (y - 2) * (up + down))
										/ ((x + y - 4) * 2);
	kernelInitInner_f<<<gridInner, block>>>(d_u, d_w, average, x, y);
	kernelInitUpDownBoundary_f<<<gridLine, block>>>(d_u, d_w, up, down, x, y);
	kernelInitLeftRightBoundary_f<<<gridColumn, block>>>(d_u, d_w, left, right, x, y);
}

//kernel of jacobi iteration process with float in iteration mode
__global__ void kernelJacobiIteration_f(int x, int y, float *d_u, float *d_w)
{
	int id = blockIdx.x * blockDim.x +threadIdx.x;
	int t = id / y;
	if ( id < x * y - 2 * y -1 && id != t * y   && id != t * y + y - 1)
	{
		id = id + y;	
		d_w[id] = (d_u[id - y] + d_u[id + y] + d_u[id - 1] 
											+ d_u[id + 1]) / 4.0;
	}	
}

//kernel of getting all epsilon between d_u and d_w
__global__ void kernelGetAllEpsilon_f(int x, int y, 
									  float *d_u, float *d_w, float *d_ep)
{
	int id = blockIdx.x * blockDim.x +threadIdx.x;
	if ( id < x * y - 2 * y -1)
	{
		id = id + y;	
		d_ep[id - y] = d_u[id] - d_w[id];
		if (d_ep[id - y] < 0) d_ep[id - y] = - d_ep[id - y];
	}	
}

//kernel of max of a num group
__global__ void kernelMaxofGroup_f(int count, float *d_ep)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	if (d_ep[id] < d_ep[id + count]) d_ep[id] = d_ep[id + count];
}

float getEpsilon_f(int x, int y, float *d_u, float *d_w)
{
	float		epsilon;
	float		*d_ep;
    cutilSafeCall(cudaMalloc((void**) &d_ep, (x * y - 2 * y) * sizeof(float)));
	int			size_block = 256;
    dim3		block(size_block);
    dim3		grid(carry(x * y - 2 * y , size_block));
	dim3		grid2(carry(x * y / 2 - y, size_block));
	kernelGetAllEpsilon_f<<<grid, block>>>(x, y, d_u, d_w, d_ep);
	//get max of d_ep
	int			count = (int)(x * y / 2 - y);
	while(count > 1)
	{
		kernelMaxofGroup_f<<<grid2, block>>>(count, d_ep);
		count = (count + 1)/ 2;
	}

	cutilSafeCall(cudaMemcpy(&epsilon, d_ep, sizeof(float), cudaMemcpyDeviceToHost));
	cutilSafeCall(cudaFree(d_ep));

	return epsilon;
}

//jacobi iteration process with float in epsilon mode
void cudaJacobiIteration_epsilon_f(int x, int y, float epsilon, 
								   float *d_u, float *d_w, int *iteration)
{
	float		*temp;	
	// setup execution parameters	
	int			size_block = 256;
    dim3		block(size_block);
    dim3		grid(carry(x * y - 2 * y , size_block));
	int			step = 0;
	float		goal_epsilon = epsilon + 1;
    
	//iteration
    while (goal_epsilon > epsilon)
	{
		step = step + 1;
		for (int i = 0; i < 1; i ++)
		{
			//execute the kernel
			kernelJacobiIteration_f<<<grid, block>>>(x, y, d_u, d_w);
			//check if kernel execution generated and error
			cutilCheckMsg("Kernel execution failed");
			temp = d_u; d_u = d_w; d_w = temp;
		}
		//epsilon
		goal_epsilon = getEpsilon_f(x, y, d_u, d_w);		
    }
	*iteration = step;
}

//jacobi iteration process with float in iteration mode
void cudaJacobiIteration_iteration_f(int x, int y, int iteration, 
									 float *d_u, float *d_w, float *epsilon)
{
	float		*temp;
	float		*d_ep;	
    cutilSafeCall(cudaMalloc((void**) &d_ep, (x * y - 2 * y) * sizeof(float)));
	// setup execution parameters	
	int			size_block = 256;
    dim3		block(size_block);
    dim3		grid(carry(x * y - 2 * y, size_block));
	dim3		grid2(carry(x * y / 2 - y, size_block));
    
	//iteration
    for (int i = 0; i < iteration; i ++)
    {
        // execute the kernel
        kernelJacobiIteration_f<<<grid, block>>>(x, y, d_u, d_w);
		// check if kernel execution generated and error
		cutilCheckMsg("Kernel execution failed");
		temp = d_u; d_u = d_w; d_w = temp;
    }
	//epsilon
	*epsilon = getEpsilon_f(x, y, d_u, d_w);	
}

//cuda jacobi with float
void cudaJacobi_f(int argc, char** argv, 
				  float left, float up, float right, float down, 
				  int x, int y,
				  float epsilon, 
				  int iteration, 
				  const char *outputfilename)
{
    if( cutCheckCmdLineFlag(argc, (const char**)argv, "device") )
        cutilDeviceInit(argc, argv);
    else
        cudaSetDevice(cutGetMaxGflopsDeviceId());

	// create and start timer
    unsigned int	timer = 0;
    cutilCheckError(cutCreateTimer(&timer));
    cutilCheckError(cutStartTimer(timer));

    // allocate host memory for matrices u
    unsigned int	size_u = x * y;
    unsigned int	mem_size_u = sizeof(float) * size_u;

    // allocate device memory
    float*			d_u;
    cutilSafeCall(cudaMalloc((void**) &d_u, mem_size_u));
    float*			d_w;
    cutilSafeCall(cudaMalloc((void**) &d_w, mem_size_u));

	//init
	cudaInitMatrix_f(x, y, left, up, right, down, d_u, d_w);

	//jacobi iteration process in epsilon mode
	if (epsilon != 0 && iteration == 0)
	{
		outputString(" - Epsilon Mode ", true);
		outputFloat(epsilon, "\tEpsilon = ",true);
		outputString("\tComputing starts...", false);

		cudaJacobiIteration_epsilon_f(x, y, epsilon, d_u, d_w, &iteration);

		outputString("Computing ends.", true);
		outputInt(iteration, "\tIteration count = ", true);
	}
	//jacobi iteration process in iteration mode
	else if (epsilon == 0 && iteration != 0)
	{
		outputString(" - Iteration Mode", true);
		outputInt(iteration, "\tIteration = ",true);
		outputString("\tComputing starts...", false);

		cudaJacobiIteration_iteration_f(x, y, iteration, d_u, d_w, &epsilon);

		outputString("Computing ends.", true);
		outputFloat(epsilon, "\tEpsilon = ", true);
	}
	else
		outputString("\tEpsilon or Iteration is wrong.", true);

    // copy result from device to host	
    float*			u = (float*) malloc(mem_size_u);
    cutilSafeCall(cudaMemcpy(u, d_u, mem_size_u,
                              cudaMemcpyDeviceToHost) );

    // stop and destroy timer
    cutilCheckError(cutStopTimer(timer));
	outputDouble(cutGetTimerValue(timer),"Processing time (ms): ", true);
    //printf("Processing time: %f (ms) \n", cutGetTimerValue(timer));
    cutilCheckError(cutDeleteTimer(timer)); 

	//output matix
	//outputMatrix_f(u, x, y, "matrix from GPU :", false);
	outputFloatMatrixtoFile(u, x, y, "matrix from GPU-CUDA with Float :", outputfilename);

    // clean up memory
    free(u);
    cutilSafeCall(cudaFree(d_u));
    cutilSafeCall(cudaFree(d_w));

    cudaThreadExit();
}

//*****************************************************************************
//Cuda jacobi with DOUBLE
//*****************************************************************************
//test double 
__global__ void 
kernelTest_d(double* u)
{
	int id = blockIdx.x * blockDim.x +threadIdx.x;
	u[id] = 0; 
}

//init non-boundarys
__global__ void 
kernelInitInner_d(double* u, double* w, double average, int x, int y)
{
	int id = blockIdx.x * blockDim.x +threadIdx.x;
	if (id < x * y)
	{
		u[id] = average; w[id] = average;
	}
}

//init up & down boundarys
__global__ void kernelInitUpDownBoundary_d(double* u, double* w, 
									       double up, double down, int x, int y)
{
	int id = blockIdx.x * blockDim.x +threadIdx.x;
		if (id < y - 1)
	{
		int upLoca = id;
		u[upLoca] = up;	w[upLoca] = up;
		int downLoca = x * y - y + 1 + upLoca;
		u[downLoca] = down;	w[downLoca] = down;
	}
}

//init left & right boundarys
__global__ void kernelInitLeftRightBoundary_d(double* u, double* w, 
											  double left, double right, int x, int y)
{
	int id = blockIdx.x * blockDim.x +threadIdx.x;
	if (id < x - 1)
	{
		int leftLoca = y + id * y;
		u[leftLoca] = left;	w[leftLoca] = left;
		int rightLoca = leftLoca - 1;
		u[rightLoca] = right; w[rightLoca] = right;
	}
}

//init with double
void cudaInitMatrix_d(int x, int y, 
				      double left, double up, double right, double down, 
				      double *d_u, double *d_w)
{
	int			size_block = 256;
	dim3		block(size_block);
	dim3		gridLine(carry(y - 1, size_block));
	dim3		gridColumn(carry(x - 1, size_block));
	dim3		gridInner(carry(x * y, size_block));

	double		average = ((x - 2) * (left + right) + (y - 2) * (up + down))
										/ ((x + y - 4) * 2);
	//outputDouble(average, "Average = ", true); getchar();
	//dim3		gridTest(1);
	//kernelTest_d<<<gridTest, block>>>(d_u);
	kernelInitInner_d<<<gridInner, block>>>(d_u, d_w, average, x, y);
	kernelInitLeftRightBoundary_d<<<gridColumn, block>>>(d_u, d_w, left, right, x, y);
	kernelInitUpDownBoundary_d<<<gridLine, block>>>(d_u, d_w, up, down, x, y);
}

//kernel of jacobi iteration process with double in iteration mode
__global__ void kernelJacobiIteration_d(int x, int y, double *d_u, double *d_w)
{
	int id = blockIdx.x * blockDim.x +threadIdx.x;
	int t = id / y;
	if ( id < x * y - 2 * y -1 && id != t * y   && id != t * y + y - 1)
	{
		id = id + y;	
		d_w[id] = (d_u[id - y] + d_u[id + y] + d_u[id - 1] 
											+ d_u[id + 1]) / 4.0;
	}	
}

//kernel of getting all epsilon between d_u and d_w
__global__ void kernelGetAllEpsilon_d(int x, int y, 
									  double *d_u, double *d_w, double *d_ep)
{
	int id = blockIdx.x * blockDim.x +threadIdx.x;
	if ( id < x * y - 2 * y )
	{
		id = id + y;	
		d_ep[id - y] = d_u[id] - d_w[id];
		if (d_ep[id - y] < 0) d_ep[id - y] = - d_ep[id - y];
	}	
}

//kernel of max of a num group
__global__ void kernelMaxofGroup_d(int count, double *d_ep)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	if (id < count)
		if (d_ep[id] < d_ep[id + count]) d_ep[id] = d_ep[id + count];
}

double getEpsilon_d(int x, int y, double *d_u, double *d_w)
{
	double		epsilon;
	double		*d_ep;
    cutilSafeCall(cudaMalloc((void**) &d_ep, (x * y - 2 * y) * sizeof(double)));
	int			size_block = 256;
    dim3		block(size_block);
    dim3		grid(carry(x * y - 2 * y , size_block));
	dim3		grid2(carry(x * y / 2 - y, size_block));
	kernelGetAllEpsilon_d<<<grid, block>>>(x, y, d_u, d_w, d_ep);
	//// allocate host memory for matrices u
 //   unsigned int	size_u = x * y - 2 * y;
 //   unsigned int	mem_size_u = sizeof(double) * size_u;
 //   double*			ep = (double*) malloc(mem_size_u);
 //   cutilSafeCall(cudaMemcpy(ep, d_ep, mem_size_u,
 //                             cudaMemcpyDeviceToHost));	
	//outputDoubleMatrixtoFile(ep, x - 2, y, "d_ep = ", "d_ep.txt");
	////outputDoubleMatrix(ep, x, y, "****d_ep = ", true);
	//free(ep);
	//get max of d_ep
	int			count = (int)(x * y / 2 - y);
	while(count > 1)
	{
		kernelMaxofGroup_d<<<grid2, block>>>(count, d_ep);
		count = (count + 1)/ 2;
	}

	cutilSafeCall(cudaMemcpy(&epsilon, d_ep, sizeof(double), cudaMemcpyDeviceToHost));
	cutilSafeCall(cudaFree(d_ep));

	return epsilon;
}

//jacobi iteration process with double in epsilon mode
void cudaJacobiIteration_epsilon_d(int x, int y, double epsilon, 
								   double *d_u, double *d_w, int *iteration)
{
	double		*temp;	
	// setup execution parameters	
	int			size_block = 256;
    dim3		block(size_block);
    dim3		grid(carry(x * y - 2 * y , size_block));
	int			step = 0;
	double		goal_epsilon = epsilon + 1;
    
	//iteration
    while (goal_epsilon > epsilon)
	{
		step = step + 1;
		
		//execute the kernel
		kernelJacobiIteration_d<<<grid, block>>>(x, y, d_u, d_w);
		//check if kernel execution generated and error
		cutilCheckMsg("Kernel execution failed");
		temp = d_u; d_u = d_w; d_w = temp;
		//epsilon
		goal_epsilon = getEpsilon_d(x, y, d_u, d_w);
		//if (step % 1000 == 0) 
		//{
		//	goal_epsilon = getEpsilon_d(x, y, d_u, d_w);
		//	outputDouble(goal_epsilon, "goal_epsilon = ", true);
		//}
    }
	*iteration = step;
}

//jacobi iteration process with double in iteration mode
void cudaJacobiIteration_iteration_d(int x, int y, int iteration, 
									 double *d_u, double *d_w, double *epsilon)
{
	double		*temp;
	double		*d_ep;	
    cutilSafeCall(cudaMalloc((void**) &d_ep, (x * y - 2 * y) * sizeof(double)));
	// setup execution parameters	
	int			size_block = 256;
    dim3		block(size_block);
    dim3		grid(carry(x * y - 2 * y, size_block));
	dim3		grid2(carry(x * y / 2 - y, size_block));
    
	//iteration
    for (int i = 0; i < iteration; i ++)
    {
        // execute the kernel
        kernelJacobiIteration_d<<<grid, block>>>(x, y, d_u, d_w);
		// check if kernel execution generated and error
		cutilCheckMsg("Kernel execution failed");
		temp = d_u; d_u = d_w; d_w = temp;
		//epsilon
		if (i % 1000 == 0) 
		{
			outputDouble(getEpsilon_d(x, y, d_u, d_w), "goal_epsilon = ", true);
		}
    }
	//epsilon
	*epsilon = getEpsilon_d(x, y, d_u, d_w);	
}

//cuda jacobi with double
void cudaJacobi_d(int argc, char** argv, 
				  double left, double up, double right, double down, 
				  int x, int y,
				  double epsilon, 
				  int iteration, 
				  const char *outputfilename)
{
	char			*s = (char*)malloc(sizeof(char) * 64);
    if (cutCheckCmdLineFlag(argc, (const char**)argv, "device"))
        cutilDeviceInit(argc, argv);
    else
        cudaSetDevice(cutGetMaxGflopsDeviceId());

	// create and start timer
    unsigned int	timer = 0;
    cutilCheckError(cutCreateTimer(&timer));
    cutilCheckError(cutStartTimer(timer));

    // allocate host memory for matrices u
    unsigned int	size_u = x * y;
    unsigned int	mem_size_u = sizeof(double) * size_u;

    // allocate device memory
    double*			d_u;
    cutilSafeCall(cudaMalloc((void**) &d_u, mem_size_u));
    double*			d_w;
    cutilSafeCall(cudaMalloc((void**) &d_w, mem_size_u));

	//init
	cudaInitMatrix_d(x, y, left, up, right, down, d_u, d_w);

	//jacobi iteration process in epsilon mode
	if (epsilon != 0 && iteration == 0)
	{
		outputString(" - Epsilon Mode ", true);
		outputDouble(epsilon, "\tEpsilon = ",true);
		outputString("\tComputing starts...", false);

		cudaJacobiIteration_epsilon_d(x, y, epsilon, d_u, d_w, &iteration);

		outputString("Computing ends.", true);
		outputInt(iteration, "\tIteration count = ", true);

		sprintf(s, "matrix from GPU-CUDA with Double - Epsilon Mode Epsilon=%lf Iteration count=%d :", 
										epsilon, iteration);
	}
	//jacobi iteration process in iteration mode
	else if (epsilon == 0 && iteration != 0)
	{
		outputString(" - Iteration Mode", true);
		outputInt(iteration, "\tIteration = ",true);
		outputString("\tComputing starts...", false);

		cudaJacobiIteration_iteration_d(x, y, iteration, d_u, d_w, &epsilon);

		outputString("Computing ends.", true);
		outputDouble(epsilon, "\tEpsilon = ", true);

		sprintf(s, "matrix from GPU-CUDA with Double - Iteration Mode Iteration=%d Epsilon=%lf :", 
										iteration, epsilon);
	}
	else
		outputString("\tEpsilon or Iteration is wrong.", true);

    // copy result from device to host	
    double*			u = (double*) malloc(mem_size_u);
    cutilSafeCall(cudaMemcpy(u, d_u, mem_size_u,
                              cudaMemcpyDeviceToHost) );

	//output matix
	//outputDoubleMatrix(u, x, y, "matrix from GPU :", false);
	outputDoubleMatrixtoFile(u, x, y, s, outputfilename);

    // stop and destroy timer
    cutilCheckError(cutStopTimer(timer));
    printf("Processing time: %f (ms) \n", cutGetTimerValue(timer));
    cutilCheckError(cutDeleteTimer(timer));    

    // clean up memory
    free(u);
    cutilSafeCall(cudaFree(d_u));
    cutilSafeCall(cudaFree(d_w));

    cudaThreadExit();
}