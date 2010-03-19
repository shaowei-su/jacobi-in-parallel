
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
		  int *n, double *epsilon, long int *step, struct boundary *b, 
		  char *outFile);



//********************************************************************************
//output
//********************************************************************************