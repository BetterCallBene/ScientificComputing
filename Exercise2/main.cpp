#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Matrix.h"
#include "CRVector.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
using namespace std;


/*
int main(int argv, char** argc)
{
	
	
	int rvecA[5] = 	{0, 0, 0,  1, 1};
	int cvecA[5] = 	{0, 1, 2,  0, 2};
	double vvecA[5] = {3, 2, 1,  1, 2};

	int rvecB[4] = 	  {0, 0, 1,  2};
	int cvecB[4] = 	  {0, 1, 1,  0};
	double vvecB[5] = {1, 2, 1,  4};

	double vvecB1[3] = {1, 0, 4};
	

	Matrix A(2, 3, 5, rvecA, cvecA, vvecA);
	Matrix B(3, 2, 4, rvecB, cvecB, vvecB);
	ColVector v(3, vvecB1);

	Matrix C = A * B;

	C.getSparse();

	Matrix D = C + C;

	D.getSparse();

	ColVector B1 = A * v;

	B1.getSparse();
	


	return 0;
}
*/

#define OToZ(i) (i-1)

double sinm(double x, double y)
{
	return sin(M_PI * y);
}

int main(int argv, char** argc)
{
	bool North, South, East, West;
	double h = 1.0/10.0; 
	int N = (int) (1.0/h + 1.0);
	int lastElm = (N-2);
	int iK = pow(lastElm, 2);

	Matrix A(iK, iK);
	Matrix Dinv(iK, iK);
	Matrix LU(iK, iK);
	
	vector<double> xv(iK);

	for(int i = 0; i < iK; i++)
		xv[i] = 0;

	ColVector f;
	ColVector x(xv);

	int curElm = 0;
	double eps = 1e-6;


	for(int i = 1; i <= lastElm; i++)
	{
		North =(i != 1);
		South = (i != lastElm);
		for(int j= 1; j <= lastElm; j++)
		{
			bool bNeedF = false;
			double K[2] = {h*j, h*i};
			West = (j != 1);
			East = (j != lastElm);

			curElm = (i-1) * lastElm + j;

			A.addValue(OToZ(curElm), OToZ(curElm), 4);
			Dinv.addValue(OToZ(curElm), OToZ(curElm), 1.0/4.0);


			if(West)
			{
				A.addValue(OToZ(curElm),  OToZ(curElm-1), -1);
				LU.addValue(OToZ(curElm),  OToZ(curElm-1), -1);
			}
			else
			{
				f.addValue(sinm(K[0], K[1]));
				bNeedF = true;
			}

			if(East)
			{
				A.addValue(OToZ(curElm),  OToZ(curElm+1), -1);
				LU.addValue(OToZ(curElm),  OToZ(curElm+1), -1);
			}
			else
			{
				f.addValue(sinm(K[0], K[1]));
				bNeedF = true;
			}

			if(South)
			{
				A.addValue(OToZ(curElm), OToZ(curElm + lastElm), -1);
				LU.addValue(OToZ(curElm), OToZ(curElm + lastElm), -1);
			}
			if(North)
			{
				A.addValue(OToZ(curElm), OToZ(curElm - lastElm), -1);
				LU.addValue(OToZ(curElm), OToZ(curElm - lastElm), -1);
			}

			if(!bNeedF)
				f.addValue(0);
			
		}
	}

	A.sort();
	LU.sort();

	double r = 0;
	do
	{
		x = Dinv * (f - LU * x);
		ColVector rsk = f - A * x;
		RowVector rsk2 = rsk.Transpose();
		r = rsk2 * rsk;

	}while(r > eps);

	printf("Residuen: %f\n\n", r);
	x.getSparse();


	return 0;
}
