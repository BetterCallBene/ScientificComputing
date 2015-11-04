#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include "Solver.h"



void example1(int neq, double t, double* y)
{
	double yin[neq];
	for(int i = 0; i < neq; i++)
		yin[i] = y[i];

	y[0] = yin[1] * yin[2];
	y[1] = (-1) * yin[0] * yin[2];
	y[2] = -0.51 * yin[0] * yin[1];
}

int main(int argv, char** argc)
{
	int neq = 3, alloc_size_yout, alloc_size_tout, timepoints = 0;
	double t0, tend, h;
	double* y0;
	double* yout;
	double* tout;

	FILE* pFile = 0;
	t0 = 0.0;
	tend = 12.0;
	h = 0.001;

	y0 = (double *) mkl_malloc(neq * sizeof(double), 64);

	y0[0] = 0;
	y0[1] = 1;
	y0[2] = 1;

	RungeKutta solver(neq, t0, tend, y0, h, &example1);
	alloc_size_yout = solver.getMemorySize();
	alloc_size_tout = solver.getMemorySizeT();
	tout = (double *) mkl_malloc(alloc_size_tout, 64);
	yout = (double *) mkl_malloc(alloc_size_yout, 64);
	
	timepoints = solver.solve(tout, yout);


	pFile = fopen("Data.out", "w");

	if ( pFile != 0)
	{

		for(int i = 0; i < timepoints; i++)
		{
			int j;
			fprintf(pFile, "%f, ", tout[i]);
			for(j = 0; j < neq-1; j++)
				fprintf(pFile, "%e, ", yout[i * neq + j]);
			
			if(i < timepoints - 1)
				fprintf(pFile, "%e\n", yout[i * neq + j]);
			else
				fprintf(pFile, "%e", yout[i * neq + j]);
		}

		fclose(pFile);
	}



	mkl_free(y0);
	mkl_free(tout);
	mkl_free(yout);

	return 0;
}