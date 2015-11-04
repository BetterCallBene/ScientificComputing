#include "CRVector.h"
#include "stdio.h"

Matrix ColVector::Transpose() const
{	
	RowVector rV(val);
	return rV;
}


Matrix RowVector::Transpose() const
{
	ColVector cV(val);
	return cV;
}

double operator*(RowVector lhs, ColVector rhs)
{
	int lhsm = lhs.size().m;
	int lhsn = lhs.size().n;
	int rhsTm = rhs.size().n;
	int rhsTn = rhs.size().m;
	if(lhsn != rhsTn)
	{
		char error_message[100] = {""};
		sprintf(error_message, "Wrong Dimension: Left (%d, %d), Right (%d, %d)", lhsm, lhsn, rhsTn, rhsTm);
		throw MatrixException(error_message);
	}
	RowVector rhsT = rhs.Transpose();
	Matrix res(lhsm, rhsTm);
	res.reserve(lhsm *rhsTm);
	Matrix::multi(lhsm, rhsTm, res, lhs, rhsT);
	return res.FirstElm();
}


ColVector operator*(Matrix m, ColVector v)
{
	int lhsm = m.size().m;
	int lhsn = m.size().n;
	int rhsTm = v.size().n;
	int rhsTn = v.size().m;
	
	if(lhsn != rhsTn)
	{
		char error_message[100] = {""};
		sprintf(error_message, "Wrong Dimension: Left (%d, %d), Right (%d, %d)", lhsm, lhsn, rhsTn, rhsTm);
		throw MatrixException(error_message);
	}

	RowVector rhsT = v.Transpose();
	ColVector res;
	res.reserve(lhsm *rhsTm);

	Matrix::multi(lhsm, rhsTm, res, m, rhsT);
	//return Matrix::add(lhs, rhs);
	return res;
}

/*

Matrix operator*(Matrix lhs, const Matrix& rhs)
{
	//return Matrix::multi(lhs, rhs);

	if(lhs.Dim.m != rhs.Dim.n)
	{
		char error_message[100] = {""};
		sprintf(error_message, "Wrong Dimension: Left (%d, %d), Right (%d, %d)", lhs.Dim.m, lhs.Dim.n, rhs.Dim.m, rhs.Dim.n);
		throw MatrixException(error_message);
	}
	Matrix rhsT = rhs.Transpose();
	int lhsm = lhs.Dim.m;
	int rhsm = rhsT.Dim.m;
	Matrix res(lhsm, rhsm);
	multi(lhsm, rhsm, res, lhs, rhsT);
	return res;
}
*/