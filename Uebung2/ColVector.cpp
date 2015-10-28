#include "ColVector.h"

ColVector::ColVector()
{

}

ColVector::~ColVector()
{

}

ColVector::ColVector(Matrix m)
{
	swap(m);
}


ColVector::ColVector(int len, double *val)
{	
	initvec(len, val);	
}

Matrix& ColVector::GetMatrix()
{
	return *this;
}

void ColVector::initvec(int len, double* val)
{
	this->Dim.m = len;
	this->Dim.n = 1;
	this->row.reserve(len);
	this->col.reserve(len);
	for(int i = 0; i < len; i++)
	{
		this->row.push_back(i);
		this->col.push_back(0);
	}

	this->val.assign(val, val + len);
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


ColVector operator+(ColVector lhs, const ColVector& rhs)
{	
	Matrix vlhs = lhs;
	Matrix vrhs = rhs;
	return Matrix::add(vlhs, vrhs);
}

ColVector operator*(ColVector lhs, const ColVector& rhs)
{
	Matrix vlhs = lhs;
	Matrix vrhs = rhs;
	 

	return Matrix::multi(vlhs, vrhs);
}