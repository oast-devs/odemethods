#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#define VECTORIZATION_32BYTE
#ifdef VECTORIZATION_32BYTE
	#define ALIGN_SIZE 4
#else
	#define ALIGN_SIZE 3
#endif

#define nextindex(i) ALIGN_SIZE * (i)

int ODE1(double x0, double t, double step, double **output, double(*fnc_ptr)(double, double));

int ODE1_3D(double x0[3], double t, double stepsize, double **output, void(*fnc_ptr)(double, double *, double *));
