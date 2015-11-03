#ifndef __MATRIX__
#define __MATRIX__

#include <vector>
#include <exception>

struct _dim{
	int m;
	int n;
};

struct MatrixException : public std::exception
{
private:
	char * error_message; 
public:
	MatrixException(char* error_message);

	const char * what () const throw ()
  	{
    	return error_message;
  	}
};


class Matrix
{
	//void internal_transpose();
public:
	friend Matrix operator*(Matrix& other1, Matrix& other2);
	friend Matrix operator+(Matrix& other1, Matrix& other2);
	friend Matrix operator-(Matrix& other1, Matrix other2);
	friend Matrix operator*(Matrix lhs, double alpha);
	friend Matrix operator*(double alpha, Matrix lhs);
protected:
	static void multi(int lhsm, int rhsm, Matrix& res, Matrix& lhs, Matrix& rhsT, bool nosparse = true);
	static void add(int m, int n, Matrix& lhs, Matrix& rhs, Matrix& res, bool nosparse = true);
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
	void shrink_to_fit();
protected:
	virtual std::vector<int> getSubCol(int start, int len);
	std::vector<double>getSubVal(int start, int len);
	void swap(const Matrix& other);
	void multiSkalar(int alpha);
public:
	virtual void sort();
	void assign(int m, int n, int nnz, int*row, int* col, double *val);
	virtual struct _dim size() const;
	virtual void getcsr(int* dsp, int* cnt) const;
	virtual void addValue(int irow, int icol, double dval);
	int getNNZ() const;
	int getCapacity() const;
	virtual void getSparse();
	void reserve(int nnz);
	double FirstElm()
	{
		return val[0];
	}
//	void Transpose();
private:
	static void bubblesort(int first, int last, std::vector<int>* row, 
		std::vector<int>* col, std::vector<double>* val);
public:
	//static void Transpose(const Matrix& source, Matrix& dist);
	virtual Matrix Transpose() const;
public:
	Matrix();
	Matrix(int m, int n);
	Matrix(int m, int n, int nnz, int* row, int* col, double* val);
	~Matrix();
};

#endif