#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <string>
#include <math.h>
#include <float.h>
#include <assert.h>

using namespace std;
const int  MAXITER = 100, NX = 5, NY = 5;
const int m = NX * NY - 2 * (NX + NY) + 4, BSR_NNZ = m * m;
const double w = 0.9, dx = 1. / (NX - 1), dy = 1. / (NY - 1), eps = 1.e-05, zero = 0.0, one = 1.0, half = 0.5;
void init(double[][m], double[][m], double[][m], double[][m], double[]);
void ShowMatrix(double[][m], string);
void ShowVector(double[], string);
void InitmatrixF(double[]);
void Input_rights(double[]);
double MAXVAL(double[], double[]);
void REL(double[], double[], double[][m], double[][m], double[][m], double[]);

int main()
{
	int i, j, k = 0;
	double M[m][m] = { NAN }, d[m] = { NAN }, L[m][m], D[m][m], R[m][m], x[m] = { NAN }, xtmp[m] = { NAN }, criteria = zero;
	init(M, L, R, D, d);
	ShowMatrix(M, "Matrix of coefficents:");
	//ShowVector(d, "Vector pravoy chasti");
	//ShowMatrix(L, "Lower triangle matrix");
	//ShowMatrix(R, "Upper triangle matrix");
	//ShowMatrix(D, "Diag matrix");
	for (i = 0; i < m; i++)
	{
		x[i] = zero;
		xtmp[i] = zero;
	}
	cout << endl << "Solution of the SLAE by Successive Over-relaxation method with eps = " << setprecision(5) << eps << endl;
	do
	{
		REL(x, xtmp, L, R, D, d);
		criteria = MAXVAL(x, xtmp);
		for (i = 0; i < m; i++)
			xtmp[i] = x[i];
		k += 1;
		cout << scientific << setprecision(14) << "Solution(" << k << "), Critetia = " << criteria << endl;
	} while ((criteria > eps || criteria == zero) && k <= MAXITER);
	if (k > MAXITER)
		cout << endl << "Failed to convergence - maximum iterations was exceeded.." << endl;
	else
	{
		cout <<endl << "DONE !!!" << endl << endl;
		cout << "Exit criteria is:" << criteria << endl << endl;
		cout << "Relaxation parameter:" << fixed << setprecision(2) << " " << w << endl << endl;
		cout << "Solution vector:" << endl;
		for (i = 0; i < m; i++)
			cout << "x(" <<  i + 1 <<") = " << fixed << setprecision(4) << x[i] << endl;
	}
	cin.get();
	return EXIT_SUCCESS;
}
void init(double M[][m], double  L[][m], double  R[][m], double  D[][m], double d[])
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			M[i][j] = zero;
			L[i][j] = zero;
			R[i][j] = zero;
			D[i][j] = zero;
		}
		d[j] = zero;
	}
	double A[BSR_NNZ];
	InitmatrixF(A);
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
			M[i][j] = A[(i * (m)+j)];
	}
	Input_rights(d);
	ShowVector(d, "Vector pravoy chasti");
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < m; j++)
		{
			if (i < j)
				R[i][j] = M[i][j] / M[i][i];
			if (i > j)
				L[i][j] = M[i][j] / M[i][i];
			if (i == j)
				D[i][j] = M[i][j] / M[i][j];
		}
		d[i] = d[i] / M[i][i];
	}
}
void ShowMatrix(double M[][m], string t)
{
	int i, j;
	cout << t << endl;
	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < m; ++j)
		cout << setw(7) << fixed << setprecision(2) << M[i][j];
		cout << endl;
	}
	cout << endl ;
}
void ShowVector(double d[], string t)
{
	int i;
	cout << t << endl;
	for (i = 0; i < m; i++)
	cout << setw(7) << fixed << setprecision(2) << d[i];
	cout << endl << endl;
}
double MAXVAL(double v1[], double v2[])
{
	int i;
	double maxval = zero;
	for (i = 0; i < m; i++)
		if (fabs(v2[i] - v1[i]) > maxval)
			maxval = fabs(v2[i] - v1[i]);
	return fabs(maxval);
}
void REL(double x[], double xt[], double L[][m], double R[][m], double D[][m], double d[])
{
	int i, j;
	double sumL, sumR;

	for (i = 0; i < m; i++)
	{
		sumL = zero; sumR = zero;
		for (j = 0; j < m; j++)
		{
			if (j > i)
				sumR = sumR + R[i][j] * xt[j];
			if (j < i)
				sumL = sumL + L[i][j] * x[j];
		}
		x[i] = (one - w) * xt[i] + w * (d[i] - sumL - sumR);
	}
	cout << endl;
}
void InitmatrixF(double block1[])
{
	for (int i = 0; i < BSR_NNZ; i++)
	{
		block1[i] = 0.;
	}
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (i == j)
			{
				block1[i * (m)+j] = -4.0;
				if (i * (m)+j != BSR_NNZ - 1)
				{
					if ((i + 1) % (NX - 2) != 0)
					{
						block1[(i * (m)+j) + 1] = 1.;
					}
					if ((m - j) > (NX - 2))
					{
						block1[(i * (m)+j) + (NX - 2)] = 1.;
					}
				}
				if (i * (m)+j != 0)
				{
					if ((i + (NX - 2)) % (NY - 2) != 0)
					{
						block1[(i * (m)+j) - 1] = 1.;
					}
					if (j > (NX - 3))
					{
						block1[(i * (m)+j) - (NY - 2)] = 1.;
					}
				}
			}
		}
	}
}
void Input_rights(double rights[m])
{
	for (int i = 0; i < m; i++) rights[i] = 0.;
	for (int i = 2; i <= NX - 1; i++)
	{
		for (int j = 2; j <= NY - 1; j++)
		{
			if (i - 1 == 1)
			{
				rights[(i - 2) * (NX - 2) + (j - 2)] -= -(dx * (j - 1)) * (dy * (j - 1));
			}
			if (j - 1 == 1)
			{
				rights[(i - 2) * (NX - 2) + (j - 2)] -= (dx * (i - 1)) * (dx * (i - 1));
			}
			if (i + 1 == NX)
			{
				rights[(i - 2) * (NX - 2) + (j - 2)] -= 1. - (dy * (j - 1)) * (dy * (j - 1));
			}
			if (j + 1 == NY)
			{
				rights[(i - 2) * (NX - 2) + (j - 2)] -= (dy * (i - 1)) * (dx * (i - 1)) - 1.;
			}
		}
	}
}