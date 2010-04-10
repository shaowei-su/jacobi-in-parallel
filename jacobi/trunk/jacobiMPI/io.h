
//********************************************************************************
//input
//********************************************************************************

//********************************************************************************
//get input arguments. there is 3 ways -
//	way 1: read arguments from default .txt file named "input.txt";(no arguments)
//	way 2: read arguments from certain .txt file whose name is the first arguments 
//			after command;(1 argument)
//	way 3: read atguments directly from arguments.(9 arguments)
int input(int argc, char* argv[], 
		  int *n, double *epsilon, long *step, struct boundary *b,
		  char *outFile);



//********************************************************************************
//output
//********************************************************************************

//get output directory name and creat the directory
char* getOutDir(int n, double epsilon, struct boundary b, 
				long step, char *outFile);


void outMatrix1DtoF(const double *m, const int n, const char *dir);


void outMatrix2DtoF(double **m, const int n, const char *dir);


void outLog(int n, double epsilon, 
			long step, struct boundary b, 
			double nTime1, double nTime2, double nTime3, 
			char *outFile, char *dir);