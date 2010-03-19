// stdafx.h : 标准系统包含文件的包含文件，
// 或是经常使用但不常更改的
// 特定于项目的包含文件
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>



// TODO: 在此处引用程序需要的其他头文件

//********************************************************************************
//user-added .h file
#include <math.h>
#include <malloc.h>
#include <windows.h>
#include <time.h>
#include <process.h>

//********************************************************************************
//inside-project .h file
#include <common.h>
#include <io.h>
#include <jacobiSerial_1D.h>
#include <jacobiSerial_2D.h>

//********************************************************************************
//constant variable
#define DEFAULT_INPUT_FILE	"input.txt"		//default input .txt filename
#define MUL					100000			//time accuracy control numer


//********************************************************************************
//structure define
struct boundary
{
	double		left;
	double		up;
	double		right;
	double		down;
};
