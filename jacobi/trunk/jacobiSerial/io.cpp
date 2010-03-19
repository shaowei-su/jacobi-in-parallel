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
		  int *n, double *epsilon, long *step, struct boundary *b, 
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
			printf("Cannot open given input file - %s.\n", argv[1]);
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
		*step = atoi(argv[3]);
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

//********************************************************************************
//get output directory name and creat the directory
char* getOutDir(int n, double epsilon, struct boundary b, 
				long step, char *outFile)
{	
	char *s = (char *)malloc(sizeof(char) * 128);
	sprintf(s, "%s_n%de%.2lfs%dlurd%.2lf%.2lf%.2lf%.2lf", 
		outFile, n, epsilon, step, b.left, b.up, b.right, b.down);
	//mkdir(s);
	return s;
}

//********************************************************************************
//
void outMatrix1DtoF(const double *m, const int n, const char *dir)
{
	char *filename = (char *)malloc(sizeof(char) * 256);
	sprintf(filename, "%s\\%s_m.txt", dir, dir);
	
	FILE *fp;
	if((fp = fopen(filename, "at")) == NULL)
	{
        printf("cannot open %s \n", filename);
		return;
	}
	else
	{
		for(int i = 0; i < n; i++)		
		{
			for(int j = 0; j < n; j++)			
				fprintf(fp, "%lf ", m[i * n + j]);
			fprintf(fp, "\n");
		}
		fclose(fp);
	}
	return;
}

//********************************************************************************
//
void outLog(int n, double epsilon, 
			long step, struct boundary b, 
			double nTime1, double nTime2, double nTime3, 
			char *outFile, char *dir)
{
	char *filename = (char *)malloc(sizeof(char) * 256);
	sprintf(filename, "%s\\%s_log.txt", dir, dir);

	time_t rawtime;   
	struct tm * timeinfo;   	  
	time(&rawtime);   
	timeinfo = localtime ( &rawtime );

	FILE *fp;
	if((fp = fopen(filename, "at")) == NULL)
	{
        printf("cannot open %s \n", filename);
		return;
	}
	else
	{
		fprintf(fp, "Jacobi Serial - %s", outFile);
		fprintf(fp, "--%s", asctime (timeinfo));
		fprintf(fp, "--N = %d, Epsilon = %lf, Step = %ld\n", n, epsilon, step);
		fprintf(fp, "--Boundary - left = %lf, right = %lf, up = %lf, down = %lf\n", 
			b.left, b.right, b.up, b.down);
		fprintf(fp, "--(Time/s)Init=%lf, Computing=%lf, Data-saving=%lf, Total=%lf\n", 
			nTime1, nTime2, nTime3, nTime1 + nTime2 + nTime3);
		fprintf(fp, "data file - %s_m.txt\n\n", dir);
		fclose(fp);
	}
	return;
}
