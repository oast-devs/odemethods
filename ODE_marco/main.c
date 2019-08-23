#include <stdio.h>
#include <stdlib.h>
#include "ODE_methods.h"

double fnc_to_be_integrated(double t, double x)
{
	return x;
}

int main() 
{
	double * output;

	int nsteps = ODE1(0.001, 10, 0.0001, &output, &(fnc_to_be_integrated));

	if(nsteps < 0) {
		fprintf(stderr, "something went horribly wrong\n");
		exit(EXIT_FAILURE);
	}
	
	fprintf(stdout, "Integration completed. Here follows the integration vector:\n");
	for(int i = 0; i < nsteps; i++) {
		fprintf(stdout, "%.4f ", output[i]);
	}
	fprintf(stdout, "\n");
	
	exit(EXIT_SUCCESS);

}
