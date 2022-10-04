#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <iomanip>
#include <string>
#include <math.h>
#include <float.h>

using namespace std;
const int p = 8, msz = 10;

void nanofication(double Matrix[][msz])
{
	int i, j;
	for (i = 0; i < msz; i++)
	{
		for (j = 0; j < msz; j++) Matrix[i][j] = NAN;
	}
}

void read_from_the_disk(double A[][msz], int msz)
{
	int i, j;
	ifstream in("GAUSS_DATA_22");
	if (!in)
	{
		cerr << "Cannot open file" << endl;
	}
	for (i = 0; i < msz; i++) {
		for (j = 0; j < msz; j++) {
			in >> A[i][j];
		}
	}
	cout << "Finished...." << endl;
	in.close();
}

double find_norm(double A[][msz], int msz)
{
	double n = 0;
	for (int i = 0; i < msz; i++)
	{
		for (int j = 0; j < msz; j++)
		{
			n += A[i][j] * A[i][j];
		}
	}
	double norm = sqrt(n);
	return norm;
}

void ShowMatrix(double Matrix[][msz], int msz, int nsz, string str)
{
	cout << endl << str << endl;
	for (int i = 0; i < msz; i++)
	{
		for (int j = 0; j < nsz; j++)
		{
			cout << " " << setw(p) << fixed << setprecision(3) << Matrix[i][j];
		}
		cout << endl;
		cout << endl;
	}
}
void create_matrix1(double E[][msz], int msz)
{
	for (int i = 0; i < msz; i++)
	{
		for (int j = 0; j < msz; j++)
		{
			E[i][j] = 0;
		}
	}
	E[0][0] = 1;
	for (int i = 0; i < msz; i++)
	{
		for (int j = 0; j < msz; j++)
		{
			if (i = j) E[i][j] = 1;
		}
	}
}
void ForwardSubstitution(double A[][msz], double E[][msz], int N)
{
	for (int k = 0; k < N; k++)
	{
		double temp;
		temp = A[k][k];
		for (int j = 0; j < N; j++)
		{
			A[k][j] /= temp;
			E[k][j] /= temp;
		}

		for (int i = k + 1; i < N; i++)
		{
			temp = A[i][k];
			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}
}
void BackSubstitution(double A[][msz], double E[][msz], int N)
{
	double temp;
	for (int k = N - 1; k > 0; k--)
	{
		for (int i = k - 1; i >= 0; i--)
		{
			temp = A[i][k];
			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}
}
void check(double A[][msz], double inverseA[][msz], double C[][msz]) {
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			C[i][j] = 0;
			for (int k = 0; k < 10; k++)
				C[i][j] += A[i][k] * inverseA[k][j];
		}
	}
}

void Saving_A_inside_B(double A[10][10], double B[10][10], int msz, int nsz)
{
	for (int i = 0; i < msz; i++)
	{
		for (int j = 0; j < nsz; j++)
		{
			B[i][j] = A[i][j];
		}
	}
}

int main()
{
	double A[msz][msz];
	double E[msz][msz];
	double C[msz][msz];
	double B[msz][msz];
	nanofication(A);
	nanofication(E);
	nanofication(C);
	nanofication(B);
	double fnorm_of_A, fnorm_of_Ax;
	read_from_the_disk(A, msz);
	ShowMatrix(A, msz, msz, "Matrix A:");
	Saving_A_inside_B(A, B, msz, msz);
	fnorm_of_A = find_norm(A, msz);
	create_matrix1(E, msz);
	ShowMatrix(E, msz, msz, "Matrix E just created:");
	ForwardSubstitution(A, E, msz);
	ShowMatrix(A, msz, msz, "Matrix A after forwardsubstitution:");
	ShowMatrix(E, msz, msz, "Matrix E after forwardsubstitution:");
	BackSubstitution(A, E, msz);
	ShowMatrix(A, msz, msz, "Matrix A after backsubstitution:");
	ShowMatrix(E, msz, msz, "Matrix E after backsubstitution (A inversed):");
	fnorm_of_Ax = find_norm(E, msz);
	cout << endl;
	cout << "Matrix A norm =" << " " << fnorm_of_A << endl;
	cout << "Matrix Ax norm =" << " " << fnorm_of_Ax << endl;
	cout << "Chislo obuslovlennosti =" << " " << fnorm_of_A * fnorm_of_Ax;
	cout << endl;
	check(B, E, C);
	ShowMatrix(C, msz, msz, "Matrix Check");
	cin.get();
	return 0;
}
