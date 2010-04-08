

//*****************************************************************************
//Serial jacobi with FLOAT
//*****************************************************************************
void serialJacobi_f(float left, float up, float right, float down, 
					int x, int y, 
					float epsilon, 
					int iteration, 
					const char *outputfilename);


//*****************************************************************************
//Serial jacobi with DOUBLE
//*****************************************************************************
void serialJacobi_d(double left, double up, double right, double down, 
					int x, int y, 
					double epsilon, 
					int iteration, 
					const char *outputfilename);