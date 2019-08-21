#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

int ODE1(double x0, double t, double step, double **output, double(*fnc_ptr)(double, double));
