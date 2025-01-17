#ifndef __ROWVECTOR__
#define __ROWVECTOR__

#include "Matrix.h"

class RowVector : public Matrix
{
public:
	RowVector& operator=(const RowVector& other);
	RowVector& operator=(const Matrix& other);
	friend RowVector operator*(RowVector other1, const RowVector& other2);
	friend RowVector operator+(RowVector other1, const RowVector& other2);
	friend RowVector operator*(RowVector lhs, double alpha);
	friend RowVector operator*(double alpha, RowVector lhs);

public:
	RowVector();
	~RowVector();
	RowVector(int len, double *val);
	RowVector(Matrix m);
protected:
	std::vector<int> getSubCol(int start, int len);
	void swap(const Matrix& other);
public:
	void getcsr(int* dsp, int* cnt) const;
	//ColVector Transpose();
};

#endif