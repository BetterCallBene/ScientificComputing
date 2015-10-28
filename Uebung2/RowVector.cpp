#include "RowVector.h"

RowVector::RowVector()
{

}

RowVector::~RowVector()
{

}

RowVector::RowVector(Matrix m)
{
	swap(m);
}


RowVector::RowVector(int len, double *val)
{	
	initvec(len, val);	
}

void RowVector::initvec(int len, double* val)
{
	this->Dim.n = len;
	this->Dim.m = 1;
	this->row.reserve(len);
	this->col.reserve(len);
	for(int i = 0; i < len; i++)
	{
		this->col.push_back(i);
		this->row.push_back(0);
	}

	this->val.assign(val, val + len);
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


RowVector operator+(RowVector lhs, const RowVector& rhs)
{	
	Matrix vlhs = lhs;
	Matrix vrhs = rhs;
	return Matrix::add(vlhs, vrhs);
}

RowVector operator*(RowVector lhs, const RowVector& rhs)
{
	Matrix vlhs = lhs;
	Matrix vrhs = rhs;
	 

	return Matrix::multi(vlhs, vrhs);
}