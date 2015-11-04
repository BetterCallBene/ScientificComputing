#include <stdlib.h>
#include <stdio.h>

#define MIN(X,Y) (X < Y ? X : Y)
#define MAX(X,Y) (X > Y ? X : Y)


/*
int main(int argv, char** argc)
{
	int len = 10;
	double res[len];
	int rescol[len];

	int colmin[3] = {1, 4, 8};
	double valmin[3] = {5, 3, 1};

	int colmax[6] = {0, 2, 3, 4, 5, 8};
	double valmax[6] = {7, 4, 5, 6, 7, 10};

	for(int i = 0; i < len; i++)
	{
		res[i] = 0;
		rescol[i] = 0;
	}

//	skalar_product(3, colmin, valmin, 8, colmax, valmax, res);

//	printf("%f\n", res);

	addition(len, 3, colmin, valmin, 6, colmax, valmax, rescol, res);


	return 0;
}
*/

double frand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int addition(int len, int nnz1, int* col1, double* val1, 
	int nnz2, int* col2, double* val2, int* colres, double* valres )
{
	int* mincol;
	int* maxcol;
	double* minval;
	double* maxval; 
	int lastmaxElm = 0;
	int lastminElm = 0;
	int ind = 0;
	int minnnz = MIN(nnz1, nnz2);
	int maxnnz = MAX(nnz1, nnz2);

	int indmax = 0;
	int indmin = 0;

	if(nnz1 == 0 && nnz2 == 0)
		return 0;

	if(minnnz == nnz1)
	{
		mincol = col1;
		minval = val1;

		maxcol = col2;
		maxval = val2;

	}
	else
	{
		mincol = col2;
		minval = val2;

		maxcol = col1;
		maxval = val1;
	}

	
	for(int i = 0; i < len; i++)
	{
		bool hasMinValue = false;
		bool hasMaxValue = false;
		int tmpMinCol = len;
		double tmpMinVal =0;
		int tmpMaxCol = len;
		double tmpMaxVal =0;

		if (indmin >= minnnz && indmax >= maxnnz)
			break;
		if(indmin < minnnz)
		{
			tmpMinCol = mincol[indmin];
			if(tmpMinCol == i)
			{
				hasMinValue = true;
				tmpMinVal = minval[indmin];
				indmin = indmin +1;
			}
		}
		if(indmax < maxnnz)
		{
			tmpMaxCol = maxcol[indmax];
			if(tmpMaxCol == i)
			{
				hasMaxValue = true;
				tmpMaxVal = maxval[indmax];
				indmax = indmax  +1;
			}
		}
		if (!(hasMinValue || hasMaxValue))
			continue;

		if(hasMinValue && hasMaxValue)
		{
			colres[ind] = tmpMinCol;
			valres[ind] = tmpMaxVal + tmpMinVal;
			ind = ind + 1;
		}
		else if(hasMaxValue)
		{
			colres[ind] = tmpMaxCol;
			valres[ind] = tmpMaxVal;
			ind = ind + 1;
		}
		else if(hasMinValue)
		{
			colres[ind] = tmpMinCol;
			valres[ind] = tmpMinVal;
			ind = ind + 1;
		}

	}

	return ind;
}


double skalar_product(int nnz1, int* col1, double* val1, 
	int nnz2, int* col2, double* val2)
{
	int lastmax = 0;
	int * mincol;
	double * minval;
	int * maxcol;
	double * maxval;

	double tmpVal = 0; 
	
	int minnnz = MIN(nnz1, nnz2);
	int maxnnz = MAX(nnz1, nnz2);

	if(minnnz == nnz1)
	{
		mincol = col1;
		minval = val1;

		maxcol = col2;
		maxval = val2;

	}
	else
	{
		mincol = col2;
		minval = val2;

		maxcol = col1;
		maxval = val1;
	}

	for(int i = 0; i < minnnz; i++)
	{
		bool found = false;
		int mincoli = mincol[i];
		int maxcolj = 0;
		int maxvalj = 0;
/*		Vorsichtig zero oder one based */
		int j;
		int jmax = MIN(mincoli, maxnnz);
		for(j= lastmax; j <= jmax; j++)
		{
			maxcolj = maxcol[j];
			if(maxcolj > mincoli)
				break;
			else if(maxcolj == mincoli)
			{
				lastmax = j;
				found = true;
				break;
			}
		}

		if (found)
		{
			tmpVal = tmpVal + maxval[j] * minval[i];
		}
	}

	return tmpVal;
}
