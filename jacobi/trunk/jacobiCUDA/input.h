

//input parameters from a .txt file named 'filename'
bool input_f(const char *filename,
			 float *left, float *up, float *right, float *down,
			 int *x, int *y, 
			 float *epsilon, 
			 int *iteration, 
			 char *outputfilename);


//input parameters(double) from a .txt file named 'filename'
bool input_d(const char *filename, 
			 double *left, double *up, double *right, double *down,
			 int *x, int *y, 
			 double *epsilon, 
			 int *iteration, 
			 char *outputfilename);
