//common method

#include "stdafx.h"

//********************************************************************************
double getCostTime(LARGE_INTEGER start, LARGE_INTEGER end)
{
	LARGE_INTEGER nFrequency;
	QueryPerformanceFrequency(&nFrequency);
	double time = (double)MUL * (end.QuadPart - start.QuadPart) 
						/ nFrequency.QuadPart;
	return time / MUL;

}