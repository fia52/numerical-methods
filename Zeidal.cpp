#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <string>
#include <math.h>
#include <float.h>
#include <algorithm>
using namespace std;

const int n = 10, p = 10, MAXITER = 100;
const double eps = 1.e-03, one = 1.0, zero = 0.0, half = 0.5;

void matrix_read(double[][n], double[], int, int);
void ShowVector(double[], int, string);
void ShowMatrix(double[][n], int, int, string);
double MAXVAL(double[], double[], int);
void Zeidel(double[][n], double[], double[], double[], int);
double norm(double[][n], int);
bool check(double[][n], int);

int main()
{
	int msz = n, nsz = n, i, j, k = 0;
	double A[n][n], B[n][n], c[n], x[n], xtmp[n], vb[n];
	double eps2, criteria = one, tmp, normb;
	for (i = 0; i < n; i++) 
	{
		c[i] = zero;
		x[i] = zero;
		xtmp[i] = zero;
		for (j = 0; j < n; j++)B[i][j] = zero;
	}
	matrix_read(A, vb, msz, nsz);
	cout << "THE Seidel METHOD" << endl << endl << endl;
	cout << "Matrix and Vector data:" << endl;
	cout << endl;
	ShowMatrix(A, msz, nsz, "Matrix A data : ");
	cout << endl;
	ShowVector(vb, nsz, "Vector vb data : ");
	cout << endl;
	cout << endl;
	for (i = 0; i < msz; i++)
	{
		c[i] = vb[i] / A[i][i];
		for (j = 0; j < nsz; j++) B[i][j] = -A[i][j] / A[i][i];
		B[i][i] = zero;
	}
	normb = norm(B, msz);
	cout << "matrix B inf norm =" << normb << endl;
	//
	if (!check(A, nsz)) cout << endl << "Matrix [A]: a full diagonal predominance wasn't found!" << endl;
	else cout << endl << "Matrix [A]: full diagonal predominance" << endl;
    //
	if (normb > one || normb <= half) eps2 = eps;
	else eps2 = eps * (one - normb) / normb;
	//
	cout << endl << "SOLUTION OF SLAE using the Zeidel metod with eps = " << eps2 << endl << endl;
	do {
		Zeidel(B, x, c, xtmp, nsz);
		criteria = fabs(MAXVAL(x, xtmp, nsz));
		for (i = 0; i < nsz; i++) { x[i] = xtmp[i]; }
		k += 1;
		cout << setw(p) << scientific << setprecision(14) << "SOLUTION (" << k << "), CRITERIA = " << criteria << endl;
	} while ((criteria > eps2 || criteria == zero) && k <= MAXITER);
	//
	if (k > MAXITER) cout << "FAILED to convergence: maximum iterations was exceeded!" << endl;
	else
	{
		cout << endl << "FINISHED!" << endl;
		cout << "iterations was made= " << k << endl;
		cout << "exit criteria = " << criteria << endl;
		cout << "solution is : " << endl << endl;
		for (i = 0; i < nsz; i++)cout << "x(" << setw(2) << i + 1 << ") =" << fixed << setprecision(3) << setw(p) << x[i] << endl;
		//
		cin.get();
		return EXIT_SUCCESS;
	}
}

void matrix_read(double matrix[][n], double vector[], int msz, int nsz)
{
	int i, j;
	ifstream is("GAUSS_DATA_22");

	if (!is) cerr << "Cannot open file" << endl;

	for (i = 0; i < msz; i++)
	{
		for (j = 0; j < nsz; j++) 
		{
			is >> matrix[i][j];
		}
	}
	for (i = 0; i < nsz; i++)is >> vector[i];
	cout << endl << "Finished reading...." << endl << endl;
	is.close();

}

void ShowMatrix(double Matrix[][n], int msz, int nsz, string str)
{
	cout << str << endl;
	for (int i = 0; i < msz; i++) {
		for (int j = 0; j < nsz; j++) {
			cout << " " << setw(p) << fixed << setprecision(4) << Matrix[i][j];
		}
		cout << endl;
	}
}

void ShowVector(double Vector[], int nsz, string str)
{
	cout << str << endl;
	for (int j = 0; j < nsz; j++) {
		cout << Vector[j] << " ";
	}
	cout << endl;
}

double norm(double B[][n], int m) {
	int i, j;
	double fnorm = zero, * tmp, sum = zero;
	tmp = new double[m];
	for (i = 0; i < m; ++i)
	{
		sum = zero;
		for (j = i; j < m; ++j)sum = sum + fabs(B[i][j]);
		tmp[i] = sum;
	}
	fnorm = *tmp, tmp + m;
	delete[] tmp;
	return fnorm;
}

bool check(double A[][n], int nsz) {
	int i, j, isdiag = true;
	double sum;
	for (i = 0; i < nsz; i++) {
		sum = zero;
		for (j = 0; j < nsz; j++) sum += fabs(A[i][j]);
		sum = sum - fabs(A[i][i]);
		if (sum >= A[i][i]) {
			isdiag = false;
			break;
		}
	}
	return isdiag;
}

void Zeidel(double mb[][n], double vx[], double vc[], double vxtmp[], int sz) 
{
	int i, j;
	double sum;
	for (i = 0; i <
		sz; i++) {
		for (j = 0; j < sz; j++) vxtmp[i] = vx[i];
	}
	for (i = 0; i < sz; i++) {
		sum = zero;
		for (j = 0; j < sz; j++) {
			sum = sum + mb[i][j] * vxtmp[j];
		}
		vxtmp[i] = sum + vc[i];
	}
}

double MAXVAL(double v1[], double v2[], int sz) {
	int i, j;
	double maxval = zero, * tmp;
	tmp = new double[sz];
	for (i = 0; i < sz; i++) tmp[i] = v1[i] - v2[i];
	maxval = *tmp, tmp + sz;
	delete[] tmp;
	return maxval;
}