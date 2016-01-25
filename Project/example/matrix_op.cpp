#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/array.hpp>
#include <cmath>

#include "drand48_generator.hpp"

using namespace std;
namespace ublas = boost::numeric::ublas;
using namespace ublas;
typedef matrix<double, row_major, array<double, 9> > Matrix;




int main(int argc, char** argv)
{
	const int n_intervals = 5;
	const int n_timepoint = n_intervals + 1;
	const int state_contr = 17;
	int n_all = n_timepoint * state_contr;
	std::vector<double> std_vec(n_all);

	srand48_generator s48;
	for(int i = 0; i < n_all; i++)
	{
		std_vec[i] = s48();
	}
	ublas::vector<double, std::vector<double>>vec(std_vec);


	cout <<"Matrixoperationen" << endl << endl; 
	array<double, 9> arrM1 = {1, 0, 0,  0, 11./(5.0*sqrt(5.0)), -2./(5.0*sqrt(5.0)), 0, 2./(5.0*sqrt(5.0)), 11./(5.0*sqrt(5.0))};
	array<double, 9> arrM2 = {3, 1, -2, 0, 11./(sqrt(5)),       11./(sqrt(5.0)),       0, -2./(sqrt(5.0)), -2./(sqrt(5.0))};

	std::vector<double> sub(arrM1.begin()+2, arrM1.begin()+6);

	Matrix A(3, 3, arrM1);
	Matrix B(3, 3, arrM2);

	cout << "Matrix A:" << A << endl;
	cout << "Matrix B:" << B << endl;

	matrix<double> res = prod(A, B);
	cout << "A * B -> " << res <<endl;

	matrix<double> subA = subrange(A, 0, 2, 0, 2);
	cout << "SubMatrix A[1..2][1..2]:" << subA <<endl;

	for (int i=0; i < n_all; i++)
	{
		cout << "vec: i=" << i  << " Value=" << vec[i] << endl;
	}

	


	for(int i =0; i < n_timepoint; i++ )
	{
		ublas::vector<double> s = 
			subrange(vec, i * state_contr + 3, i * state_contr + 7);
		cout << "subrange" << s <<endl; 
		double n2 = norm_2(s);
		for(int j = 0; j < 4; j++)
		{
			int indx = i * state_contr + 3 + j;
			vec[indx] = vec[indx] / n2;
		}
	}

	for(int i =0; i < n_timepoint; i++ )
	{
		ublas::vector<double> s = 
			subrange(vec, i * state_contr + 3, i * state_contr + 7); 
		double n2 = norm_2(s);
		cout << "Norm: " << n2 << endl;
	}

	std_vec =vec.data();
	int indx = 0;
	for(int i = 0; i < n_all; i++)
	{
		if(std_vec[i] != vec[i])
		{
			indx++;
			cout << "Error: std_vec: " << std_vec[i] << " vec: " << vec[i] << endl;
		}
	}	

	matrix<double, row_major, std::vector<double>> m = identity_matrix<double>(3);



	cout << "Identity Matrix: " << m <<endl; 

	return 0;
}