#ifndef __COLVECTOR__
#define __COLVECTOR__

#include "Matrix.h"

class ColVector : public Matrix
{
public:
	ColVector& operator=(const ColVector& other);
	ColVector& operator=(const Matrix& other);
	friend ColVector operator*(ColVector other1, const ColVector& other2);
	friend ColVector operator+(ColVector other1, const ColVector& other2);
	friend ColVector operator*(ColVector lhs, double alpha);
	friend ColVector operator*(double alpha, ColVector lhs);

public:
	ColVector();
	~ColVector();
	ColVector(int len, double *val);
	ColVector(Matrix m);
	Matrix& GetMatrix();
private:
	void initvec(int len, double* val);
};

#endif