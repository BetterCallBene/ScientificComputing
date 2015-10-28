#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include "Matrix.h"
#include "helper.h"
#include <cmath>

using namespace std;


Matrix::Matrix() 
{}

Matrix::Matrix(int m, int n)
{
	this->Dim.m = m;
	this->Dim.n = n;
}

Matrix::Matrix(int m, int n, int nnz, int* row, int* col, double* val)
{
	assign(m, n, nnz, row, col, val);
}
Matrix::~Matrix()
{
	
}

void Matrix::assign(int m, int n, int nnz, int*row, int* col, double *val)
{
	this->Dim.m = m;
	this->Dim.n = n;
	this->row.assign(row, row + nnz);
	this->col.assign(col, col + nnz);
	this->val.assign(val, val + nnz);
	this->sort();
}


void Matrix::reserve(int capacity)
{
	this->row.reserve(capacity);
	this->col.reserve(capacity);
	this->val.reserve(capacity);
}

void Matrix::shrink_to_fit()
{
	int nnz = getNNZ();
	this->row.resize(nnz);
	this->col.resize(nnz);
	this->val.resize(nnz);
}

void Matrix::swap(const Matrix& other)
{
	int nnz = other.getNNZ();
	this->Dim = other.Dim;
	this->reserve(nnz);
	this->row = other.row;
	this->col = other.col;
	this->val = other.val; 
	
}

struct _dim Matrix::size() const
{
	return Dim;
}

void Matrix::addValue(int irow, int icol, double dval)
{
	int nnz = getNNZ();
	row.push_back(irow);
	col.push_back(icol);
	val.push_back(dval);
}

int Matrix::getNNZ() const
{
	return this->row.size();
}

void Matrix::getSparse()
{
	int nnz = getNNZ();
	int capacity = getCapacity();
	printf("Sparse-Matrix -nnz: %d -dim (%d, %d) -capacity%d\n", nnz, Dim.m, Dim.n, capacity);
	for(int i = 0; i < nnz; i++)
	{
		printf("row: %d, col: %d, val:%f\n", row[i], col[i], val[i]);
	}
	printf("\n");
}

void Matrix::sort()
{
	int nnz = this->getNNZ();
	Matrix::bubblesort(0, nnz, &this->row, &this->col, &this->val);
	int first = row[0];
	int ind = 0;
	int j = 0;

	for(int i = 0; i < nnz; i++)
	{
//		printf("row: %d, col: %d\n", row[i], col[i]);
		int rowi = row[i];
		if(first < rowi)
		{
			if(ind > 1)
				Matrix::bubblesort(j, j + ind, &this->col, &this->row, &this->val);
			j = j + ind;
			ind = 1;
			first = rowi;
		}
		else
			ind++;
		
	}
	if(ind > 1)
		Matrix::bubblesort(j, j + ind, &this->col, &this->row, &this->val);

}
void Matrix::bubblesort(int first, int last, std::vector<int>* row, 
	std::vector<int>* col, std::vector<double>* val)
{
	int m = last;
	
	for (m; m>1; m--){
    	for (int i=first; i<m-1; i++){
    		int rowi  = row->at(i);
    		int rowiplus = row->at(i+1);
      		if (rowi > rowiplus){
//        		A.swap(i, i+1)
      			iter_swap(row->begin()+i, row->begin()+ (i+1));
      			iter_swap(col->begin()+i, col->begin()+ (i+1));
      			iter_swap(val->begin()+i, val->begin()+ (i+1));
      		} // ende if
    	} // ende innere for-Schleife
  	}
  
}


void Matrix::getcsr(int* dsp, int* cnt) const
{
	int nnz = getNNZ();
	int m = this->Dim.m;
	int curRow = 0;
	int ind = 0;
	int rowindex = 0;
	/*Empty*/
	if (nnz == 0)
		return;

	for(int i = 0; i < m; i++)
	{
		if(rowindex >= nnz)
		{
			cnt[i] = 0;
			dsp[i] = rowindex-1;
			continue;
		}
		cnt[i] = 0;
		dsp[i] = rowindex;
		if(i < row[rowindex])
			continue;
		
		while(rowindex <= nnz)
		{
			if(i == row[rowindex])
			{
				rowindex++;
				cnt[i]++;
			}
			else
				break;
		}

	}

}

double Matrix::norm()
{
		int nnz = getNNZ();

		double res = 0;

		for(int i = 0; i <nnz; i++)
		{
			res = res + pow(val[i], 2.0);
		}
		return sqrt(res);
}

int Matrix::getCapacity() const
{
	return row.capacity();
}
void Matrix::Transpose(const Matrix& source, Matrix& dist)
{
	//dist = source;
	int nnz = source.getNNZ();

	dist.Dim.n = source.Dim.m;
	dist.Dim.m = source.Dim.n;

	dist.reserve(nnz);

	for (int i = 0; i < nnz; i++)
	{
		//dist.row[i] = source.col[i];
		//dist.col[i] = source.row[i];
		//dist.val[i] = source.val[i];
		dist.row.push_back(source.col[i]);
		dist.col.push_back(source.row[i]);
		dist.val.push_back(source.val[i]);
	}

	dist.sort();
}

void Matrix::multiSkalar(int alpha)
{
	int nnz = getNNZ();
	for(int i = 0; i <nnz; i++)
	{
		val[i] = val[i] * alpha;
	}
}

Matrix Matrix::multi(Matrix lhs, const Matrix& rhs)
{
	Matrix rhsT;
	Matrix::Transpose(rhs, rhsT);

	int lhsm = lhs.Dim.m;
	int rhsm = rhsT.Dim.m;

	int dspLhs[lhsm];
	int cntLhs[lhsm];

	int dspRhs[rhsm];
	int cntRhs[rhsm];

	double eps = 1e-20;


	lhs.getcsr(dspLhs, cntLhs);
	rhsT.getcsr(dspRhs, cntRhs);
	Matrix res(lhsm, rhsm);
	res.reserve(lhsm * rhsm);

	int ind = 0;


	for (int i = 0; i < lhsm; i++)
	{
		int TmpdspRhs = dspRhs[i];
		int TmpcntRhs = cntRhs[i];

		if (TmpcntRhs == 0)
			continue;

		vector<int>::iterator firstColRhs 
			= rhsT.col.begin() + TmpdspRhs;
		vector<int>::iterator lastColRhs 
			= rhsT.col.begin() + TmpdspRhs + TmpcntRhs;
		vector<int> curColRhs(firstColRhs, lastColRhs);
		vector<double>::iterator firstValRhs 
			= rhsT.val.begin() + TmpdspRhs;
		vector<double>::iterator lastValRhs 
			= rhsT.val.begin() + TmpdspRhs + TmpcntRhs;
		vector<double> curValRhs(firstValRhs, lastValRhs);
		

		for(int j = 0; j < rhsm; j++)
		{
			int TmpdspLhs = dspLhs[j];
			int TmpcntLhs = cntLhs[j];

			if (TmpcntLhs == 0)
				continue;

			vector<int>::iterator firstColLhs = lhs.col.begin() + TmpdspLhs;
			vector<int>::iterator lastColLhs = lhs.col.begin() + TmpdspLhs + TmpcntLhs;
			vector<int> curRowLhs(firstColLhs, lastColLhs);

			vector<double>::iterator firstValLhs = lhs.val.begin() + TmpdspLhs;
			vector<double>::iterator lastValLhs = lhs.val.begin() + TmpdspLhs + TmpcntLhs;
			vector<double> curValLhs(firstValLhs, lastValLhs);

			double erg = skalar_product(TmpcntLhs, curRowLhs.data(), curValLhs.data(), 
				TmpcntRhs, curColRhs.data(), curValRhs.data());

			if(abs(erg) > eps)
			{
				res.addValue(j, i, erg);
			}	
		}
	}
	res.sort();
	return res;
}

Matrix Matrix::add(Matrix lhs, const Matrix& rhs)
{
	Matrix rhsTmp = rhs;
	int lhsm = lhs.Dim.m;
	int lhsn = lhs.Dim.n;
	int rhsm = rhs.Dim.m;


	int dspLhs[lhsm];
	int cntLhs[lhsm];

	int dspRhs[rhsm];
	int cntRhs[rhsm];

	lhs.getcsr(dspLhs, cntLhs);
	rhs.getcsr(dspRhs, cntRhs);

	Matrix res(lhsm, lhsn);
	res.reserve(lhsm * lhsn);

	for(int i = 0; i < lhsm; i++)
	{
		int TmpdspRhs = dspRhs[i];
		int TmpcntRhs = cntRhs[i];

		int TmpdspLhs = dspLhs[i];
		int TmpcntLhs = cntLhs[i];

		vector<int>::iterator firstColLhs = lhs.col.begin() + TmpdspLhs;
		vector<int>::iterator lastColLhs = lhs.col.begin() + TmpdspLhs + TmpcntLhs;
		vector<int> curRowLhs(firstColLhs, lastColLhs);

		vector<double>::iterator firstValLhs = lhs.val.begin() + TmpdspLhs;
		vector<double>::iterator lastValLhs = lhs.val.begin() + TmpdspLhs + TmpcntLhs;
		vector<double> curValLhs(firstValLhs, lastValLhs);

		vector<int>::iterator firstColRhs = rhsTmp.col.begin() + TmpdspRhs;
		vector<int>::iterator lastColRhs = rhsTmp.col.begin() + TmpdspRhs + TmpcntRhs;
		vector<int> curColRhs(firstColRhs, lastColRhs);

		vector<double>::iterator firstValRhs = rhsTmp.val.begin() + TmpdspRhs;
		vector<double>::iterator lastValRhs = rhsTmp.val.begin() + TmpdspRhs + TmpcntRhs;
		vector<double> curValRhs(firstValRhs, lastValRhs);

		int colRes[lhsn];
		double valRes[lhsn];

		int nnz =  addition(lhsn, TmpcntLhs, curRowLhs.data(), curValLhs.data(), 
			TmpcntRhs, curColRhs.data(), curValRhs.data(), colRes, valRes );
	/*
		addition (int len, int nnz1, int* col1, double* val1, 
	int nnz2, int* col2, double* val2, int* colres, double* valres )
	*/
		for(int j = 0; j < nnz; j++)
			res.addValue(i, colRes[j], valRes[j]);

	}
	return res;
}

Matrix& Matrix::operator=(const Matrix& other)
{
	swap(other);
	return *this;
} 

Matrix operator*(Matrix lhs, double alpha)
{
	lhs.multiSkalar(alpha);
	return lhs;
}

Matrix operator*(double alpha, Matrix lhs)
{
	lhs.multiSkalar(alpha);
	return lhs;
}


Matrix operator+(Matrix lhs, const Matrix& rhs)
{	
	return Matrix::add(lhs, rhs);
}

Matrix operator*(Matrix lhs, const Matrix& rhs)
{
	return Matrix::multi(lhs, rhs);
}


