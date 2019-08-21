#include "ODE1.h"

int ODE1(double x0, double t, double stepsize, double **output, double(*fnc_ptr)(double, double))
{
	double nsteps_tmp = t / stepsize;
	if(nsteps_tmp > (double) INT_MAX){
		fprintf(stderr, "Error: intregration step is too small\n");
		return -1;
	}
	
	int nsteps = (int) nsteps_tmp;
	
	*output = malloc(nsteps * sizeof(double));
	if(*output == NULL) {
		fprintf(stderr, "Error: system memory is unavailable\n");
		return -1;
	}
	
	(*output)[0] = x0;
	for(int i = 1; i < nsteps; i++) {
		double step = (*fnc_ptr)(i * stepsize, (*output)[i - 1]);
		(*output)[i] = (*output)[i - 1] + stepsize * step;
	}
	
	return nsteps;
}
