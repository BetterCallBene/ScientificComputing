#ifndef __MATRIX__
#define __MATRIX__

#include <vector>

struct _dim{
	int m;
	int n;
};


class Matrix
{
	//void internal_transpose();
public:
	friend Matrix operator*(Matrix other1, const Matrix& other2);
	friend Matrix operator+(Matrix other1, const Matrix& other2);
	friend Matrix operator*(Matrix lhs, double alpha);
	friend Matrix operator*(double alpha, Matrix lhs);
protected:
	static Matrix multi(Matrix other1, const Matrix& other2);
	static Matrix add(Matrix other1, const Matrix& other2);
public:
	Matrix& operator=(const Matrix& other);
protected:
	_dim Dim;
protected:
	std::vector<int> row;
	std::vector<int> col;
	std::vector<double> val;
	/*The col vector stores all column indexes of the non-zero matrix entries in row-wise order*/
//	int *col;
	/*The row vector stores all rows indexes of the non zero matrix entries*/
//	int *row;
	/*The val vector stores all values of the non-zero matrix entries in row-wise order*/
//	double *val;

private:
	void reserve(int nnz);
	void shrink_to_fit();
protected:
	void sort();
	void swap(const Matrix& other);
	void multiSkalar(int alpha);
public:
	void assign(int m, int n, int nnz, int*row, int* col, double *val);
	struct _dim size() const;
	void getcsr(int* dsp, int* cnt) const;
	void addValue(int irow, int icol, double dval);
	int getNNZ() const;
	int getCapacity() const;
	void getSparse();

	double norm();

	double FirstElm()
	{
		return val[0];
	}
//	void Transpose();
private:
	static void bubblesort(int first, int last, std::vector<int>* row, 
		std::vector<int>* col, std::vector<double>* val);
public:
	static void Transpose(const Matrix& source, Matrix& dist);
public:
	Matrix();
	Matrix(int m, int n);
	Matrix(int m, int n, int nnz, int* row, int* col, double* val);
	~Matrix();
};

#endif