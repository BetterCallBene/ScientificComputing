#include "CRVector.h"
#include <vector>
#include <stdio.h>

using namespace std;

RowVector::RowVector()
{
	this->Dim.m = 1;
}

RowVector::~RowVector()
{

}

RowVector::RowVector(Matrix m)
{
	swap(m);
}

RowVector::RowVector(vector<double> val)
{
	Dim.n = val.size();
	Dim.m = 1;

	this->val = val;
}


RowVector::RowVector(int len, double *val)
{	
	this->Dim.n = len;
	this->Dim.m = 1;
	
	this->val.assign(val, val + len);
}

void RowVector::getcsr(int* dsp, int* cnt) const
{
	dsp[0] = 0;
	cnt[0] = this->Dim.n;
}


vector<int> RowVector::getSubCol(int start, int len)
{
	vector<int> cols(len);
	int ind = 0;
	for(int i = start; i < len; i++)
		cols[ind++] = i;
	return cols;
}

void RowVector::swap(const Matrix& other)
{
	if(other.size().m != 1)
		throw MatrixException("Matrix is not a column vector");
	Matrix::swap(other);
}

void RowVector::addValue(double val)
{
	this->val.push_back(val);
	this->Dim.n = this->val.size();
}

void RowVector::getSparse()
{
	int nnz = getNNZ();
	printf("Vector: Col: %d\n", nnz);
	for(int i = 0; i < nnz; i++)
	{
		printf("%f", val[i]);
		if(i < nnz - 1)
			printf(",");
	}
	printf("\n");
}

void RowVector::addValue(int irow, int icol, double dval)
{
	addValue(dval);
}

RowVector operator-(RowVector& lhs, RowVector rhs)
{
	double alpha = -1;
	rhs.multiSkalar(alpha);
	return lhs + rhs;
}

RowVector& RowVector::operator=(const Matrix& other)
{
	swap(other);
	return *this;
}

RowVector& RowVector::operator=(const RowVector& other)
{
	swap(other);
	return *this;
} 

RowVector operator*(RowVector lhs, double alpha)
{
	lhs.multiSkalar(alpha);
	return lhs;
}

RowVector operator*(double alpha, RowVector lhs)
{
	lhs.multiSkalar(alpha);
	return lhs;
}


RowVector operator+(RowVector& lhs, RowVector rhs)
{	
	_dim DimLhs = lhs.size();
	_dim DimRhs = rhs.size();
	
	if(DimLhs.m != DimRhs.m || DimLhs.n != DimRhs.n)
	{
		char error_message[100] = {""};
		sprintf(error_message, "Wrong Dimension: Left (%d, %d), Right (%d, %d)", DimLhs.m, DimLhs.n, DimRhs.m, DimRhs.n);
		throw MatrixException(error_message);
	}

	RowVector res;
	res.reserve(DimLhs.m * DimLhs.n);

	Matrix::add(DimLhs.m, DimLhs.n, lhs, rhs, res);
	//return Matrix::add(lhs, rhs);
	return res;
}