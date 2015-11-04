#ifndef __CRVECTOR__
#define __CRVECTOR__

#include "Matrix.h"

class RowVector;

class ColVector : public Matrix
{
public:
	ColVector& operator=(const ColVector& other);
	ColVector& operator=(const Matrix& other);
	
	friend ColVector operator+(ColVector& other1, ColVector other2);
	friend ColVector operator*(ColVector lhs, double alpha);
	friend ColVector operator*(double alpha, ColVector lhs);
	friend ColVector operator*(Matrix m, ColVector v);
	friend ColVector operator-(ColVector& lhs, ColVector rhs);
public:
	ColVector();
	~ColVector();
	ColVector(std::vector<double> val);
	ColVector(int len, double *val);
	ColVector(Matrix m);
	Matrix& GetMatrix();
protected:
	std::vector<int> getSubCol(int start, int len);
	void swap(const Matrix& other);
public:
	void addValue(double val);
	void addValue(int irow, int icol, double dval);
	void getSparse();
	void sort(){};
	virtual void getcsr(int* dsp, int* cnt) const;
	virtual Matrix Transpose() const;
};


class RowVector : public Matrix
{
public:
	RowVector& operator=(const RowVector& other);
	RowVector& operator=(const Matrix& other);
	friend RowVector operator+(RowVector& other1, RowVector other2);
	friend RowVector operator*(RowVector lhs, double alpha);
	friend RowVector operator*(double alpha, RowVector lhs);
	friend double operator*(RowVector lhs, ColVector rhs);
	friend RowVector operator-(RowVector& lhs, RowVector rhs);
public:
	RowVector();
	~RowVector();
	RowVector(std::vector<double> val);
	RowVector(int len, double *val);
	RowVector(Matrix m);

protected:
	std::vector<int> getSubCol(int start, int len);
	void swap(const Matrix& other);
public:
	void addValue(double val);
	void addValue(int irow, int icol, double dval);
	void getSparse();
	void sort(){};
	virtual void getcsr(int* dsp, int* cnt) const;
	virtual Matrix Transpose() const;
};

#endif