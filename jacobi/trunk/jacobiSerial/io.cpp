//input & output

#include "stdafx.h"

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
		  char *outFile)
{
	if (argc == 1)
	{
		char		*filename = DEFAULT_INPUT_FILE;
		FILE		*fp;
		if ((fp = fopen(filename, "r")) == NULL)
		{
			printf("Cannot open defalut input file - %s.\n", filename);
			return 1;
		}
		else
		{
			fscanf(fp, "%lf %lf %lf %lf\n", 
				&(*b).left, &(*b).up, &(*b).right, &(*b).down);
			fscanf(fp, "%d\n", n);
			fscanf(fp, "%lf\n", epsilon);
			fscanf(fp, "%d\n", step);
			fscanf(fp, "%s\n", outFile);
			fclose(fp);
		}
	}
	else if (argc == 2)
	{
		FILE		*fp;
		if ((fp = fopen(argv[1], "r")) == NULL)
		{
			printf("Cannot open given input file - %s.\n", filename);
			return 2;
		}
		else
		{
			fscanf(fp, "%lf %lf %lf %lf\n", 
				&(*b).left, &(*b).up, &(*b).right, &(*b).down);
			fscanf(fp, "%d\n", n);
			fscanf(fp, "%lf\n", epsilon);
			fscanf(fp, "%d\n", step);
			fscanf(fp, "%s\n", outFile);
			fclose(fp);
		}
	}
	else if (argc == 9)
	{
		*n = atoi(argv[1]);
		*epsilon = atof(argv[2]);
		*step = stoi(argv[3]);
		(*b).left = atof(argv[4]);
		(*b).up = atof(argv[5]);
		(*b).right = atof(argv[6]);
		(*b).down = atof(argv[7]);
		strcpy(outFile, argv[8]);			
	}
	else
	{
		printf("The count of arguments is wrong. Check it. \n");
		return 3;
	}	

	return 0;
	//000000
}


//********************************************************************************
//output
//********************************************************************************