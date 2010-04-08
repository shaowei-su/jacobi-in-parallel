

//*****************************************************************************
//Cuda jacobi with FLOAT
//*****************************************************************************
void cudaJacobi_f(int argc, char** argv, 
				  float left, float up, float right, float down, 
				  int x, int y,
				  float epsilon, 
				  int iteration, 
				  const char *outputfilename);

//*****************************************************************************
//Cuda jacobi with DOUBLE
//*****************************************************************************
void cudaJacobi_d(int argc, char** argv, 
				  double left, double up, double right, double down, 
				  int x, int y,
				  double epsilon, 
				  int iteration, 
				  const char *outputfilename);