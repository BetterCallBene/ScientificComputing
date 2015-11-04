#include "Solver.h"
#include <stdio.h>
#include <cmath>
#include "mkl.h"


void add(int n, double alpha, double* x, double beta, double* y)
{
	for(int i = 0; i < n; i++)
	{
		y[i] = alpha * x[i] + beta * y[i];
	}
}

void skalar(int n, double alpha, double* y)
{
	for(int i = 0; i < n; i++)
	{
		y[i] = alpha * y[i];
	}
}

void cpy(int n, double* x, double *y)
{
	for(int i = 0; i < n; i++)
	{
		y[i] = x[i];
	}
}



void Solver::init(int neq, double t0, double tend, double* y0, double h, void(*func)(int neq, double t, double* y))
{
	this->neq = neq;
	this->t0 = t0;
	this->tend = tend;
	this->y0 = y0;
	this->h = h;
	this->func = func;
	this->timepoints = (ceil((tend - t0) / h) + 1);
	this->nelm = this->timepoints * neq;
}

int Solver::getMemorySizeT()
{
	return this->timepoints * sizeof(double);
}

int Solver::getMemorySize()
{
	return  this->nelm *  sizeof(double);
}

int Solver::solve(double *tout, double* yout)
{
	
	bool done = false;
	int neq = this->neq;
	int nelm = this->nelm;
	int ind = 0;

	double ylast[neq]; 
	double h = this->h;
	double tstep = this->t0;
	double tend = this->tend;

	double eps = 1e-10;

	double tmp, tmp1;

	for(int i = 0; i < neq; i++ )
	{
		yout[i] = this->y0[i];
		ylast[i] = this->y0[i];
	}

	for (int i = neq; i < nelm; i++)
	{
		yout[i] = 0;
	}


	if(tend <= tstep)
	{
		printf("Error");
		return 0;
	}

	tout[ind] = tstep;

	while(!done)
	{
		if(fabs(tstep + h - tend) < eps )
		{
			h = tend - tstep;
			done = true;
		}

		this->step(tstep, ylast, h);
		ind = ind + 1;
		tstep = tstep + h;

		tout[ind] = tstep;

		for (int i = 0; i < neq; i++)
			yout[ind * neq + i] = ylast[i];
		
	}
	return this->timepoints;
}

void ForwEuler::step(double tstep, double *y, double h)
{
	int neq = this->neq;
	double ylast[neq];

	for(int i = 0; i < neq; i++)
		ylast[i] = y[i];

	this->func(neq, tstep, y);

	add(neq, 1.0, ylast, h, y);

}

void RungeKutta::step(double tstep, double *y, double h)
{
	int neq = this->neq;
	
	double k1[neq];
	double k2[neq];
	double k3[neq];
	double k4[neq];

	double ylast[neq];

// k1	
//	memcpy(k1, y, neq);
	for(int i = 0; i < neq; i++)
	{
		k1[i] = y[i];
		ylast[i] = y[i];
	}
	this->func(neq, tstep, k1);
	skalar(neq, h, k1);	
//k2 
	for(int i = 0; i < neq; i++)
		k2[i] = ylast[i] + 0.5 * k1[i];  

	this->func(neq, tstep + 0.5 * h, k2);
	skalar(neq, h, k2);
//k3
	for(int i = 0; i < neq; i++)
		k3[i] = ylast[i] + 0.5 * k2[i];
	this->func(neq, tstep + 0.5 * h, k3);
	skalar(neq, h, k3);

//k4
	for(int i = 0; i < neq; i++)
		k4[i] = ylast[i] + k3[i];

	this->func(neq, tstep + h, k4);
	skalar(neq, h, k4);

	for(int i = 0; i < neq; i++)
		y[i] = ylast[i] + 1.0/6.0*(k1[i] + 2.0 * k2[i] + 2.0*k3[i] + k4[i]);

}