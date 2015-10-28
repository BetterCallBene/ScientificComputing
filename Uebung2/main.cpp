#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Matrix.h"
#include "helper.h"
#include "ColVector.h"
#include "RowVector.h"

int main(int argv, char** argc)
{
	const int dim = 8;
	const int nnz = 10;

	const double tol = 1e-5;

	double xV[dim];
	int rvecA[nnz] = 	{1, 1,  1,  2, 5, 3, 5, 5,  7,  7};
	int cvecA[nnz] = 	{2, 1,  3,  1, 2, 4, 1, 6,  7,  6};
	double vvecA[nnz] = {4, -1, -4, 1, 2, 6, 4, 10, 11, 20};
	double bV[dim] = {1, 0, 2, 2, 3, 6, 7, 0};

	Matrix A(dim, dim, nnz, rvecA, cvecA, vvecA);

	srand(time(NULL));

	for(int i = 0; i < dim; i++)
	{
		xV[i] = frand(0, 1);
		printf("%f\n", xV[i]);
	}

	ColVector cx = ColVector(dim, xV);
	ColVector cb = ColVector(dim, bV);
	Matrix x = cx.GetMatrix();
	Matrix b = cb.GetMatrix(); 


	Matrix r = b + ((-1)* A*x);
	Matrix d = r;

	do
	{
		Matrix z = A*d;
		double rskalar = (r * r).FirstElm();
		double alpha = rskalar / (d * z).FirstElm();
		x = x + alpha * d;
		double alpha_minus = (-1)* alpha;
		r = r + alpha_minus *z;

		double beta = (r * r).FirstElm() / rskalar;
		d = r + beta * d;
	}while(r.norm() < tol);
	
	return 0;
}