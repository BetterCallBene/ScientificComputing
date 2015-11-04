
#ifndef __HELPER__
#define __HELPER__
double skalar_product(int nnz1, int* col1, double* val1, 
	int nnz2, int* col2, double* val2);

int addition(int len, int nnz1, int* col1, double* val1, 
	int nnz2, int* col2, double* val2, int* colres, double* valres );

double frand(double fMin, double fMax);

#endif