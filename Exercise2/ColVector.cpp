#include <vector>
#include "CRVector.h"
#include <stdio.h>

using namespace std;

ColVector::ColVector()
{
	Dim.n = 1;
}

ColVector::~ColVector()
{

}

ColVector::ColVector(vector<double> val)
{
	Dim.m = val.size();
	Dim.n = 1;

	this->val = val;
}

ColVector::ColVector(Matrix m)
{
	swap(m);
}


ColVector::ColVector(int len, double *val)
{	
	this->Dim.m = len;
	this->Dim.n = 1;
	this->val.assign(val, val + len);	
}

Matrix& ColVector::GetMatrix()
{
	return *this;
}

void ColVector::getcsr(int* dsp, int* cnt) const
{
	for(int i = 0; i < Dim.m; i++)
	{
		dsp[i] = i;
		cnt[i] = 1;
	}
}

vector<int> ColVector::getSubCol(int start, int len)
{
	vector<int> cols(len);
	int ind = 0;
	for(int i = start; i < len; i++)
		cols[ind++] = i;
	return cols;
}


void ColVector::swap(const Matrix& other)
{
	if(other.size().n != 1)
		throw MatrixException("Matrix is not a column vector");

	Matrix::swap(other);
}

void ColVector::getSparse()
{
	int nnz = getNNZ();
	printf("Vector: Rows: %d\n", nnz);
	for(int i = 0; i < nnz; i++)
	{
		printf("%f\n", val[i]);
	}
	printf("\n");
}

void ColVector::addValue(double val)
{
	this->val.push_back(val);
	Dim.m = this->val.size();
}

void ColVector::addValue(int irow, int icol, double dval)
{
	addValue(dval);
}



ColVector& ColVector::operator=(const Matrix& other)
{
	swap(other);
	return *this;
}

ColVector& ColVector::operator=(const ColVector& other)
{
	swap(other);
	return *this;
} 

ColVector operator*(ColVector lhs, double alpha)
{
	lhs.multiSkalar(alpha);
	return lhs;
}

ColVector operator*(double alpha, ColVector lhs)
{
	lhs.multiSkalar(alpha);
	return lhs;
}

ColVector operator-(ColVector& lhs, ColVector rhs)
{
	double alpha = -1;
	rhs.multiSkalar(alpha);
	return lhs + rhs;
}


ColVector operator+(ColVector& lhs, ColVector rhs)
{	
	_dim DimLhs = lhs.size();
	_dim DimRhs = rhs.size();
	
	if(DimLhs.m != DimRhs.m || DimLhs.n != DimRhs.n)
	{
		char error_message[100] = {""};
		sprintf(error_message, "Wrong Dimension: Left (%d, %d), Right (%d, %d)", DimLhs.m, DimLhs.n, DimRhs.m, DimRhs.n);
		throw MatrixException(error_message);
	}

	ColVector res;
	res.reserve(DimLhs.m * DimLhs.n);

	Matrix::add(DimLhs.m, DimLhs.n, lhs, rhs, res);
	return res;
}