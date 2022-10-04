// Программа находит решение дифференциального уравнения вида y'' + p(x)y' + q(x)y = f(x) а именно y'' + xy' - 2y = -1.5
#include <iostream>
#include <stdlib.h>
#include <iomanip>

using namespace std;
const int N = 100;
double a = 0; // Начало отрезка. 
double b = 1; // Конец отрезка.
double Y0 = 2; // Значение в начале.
double Yn = 2 * sqrt(2); // Значение на конце.
double q = -2; // Значение q(x) - константы.
double f = -1.5; // Значение f(x) - константы.
//
void Thomas(double C[][N], double[], double[], int);
void zeroing(double Matrix[][N]);
void ShowMatrix(double A[][N], double[], double[], int);
void Proverka(double matrix[][N], double[N], int);
void Input_rights(double rights[], int, double, double[], double, double, double, int);
void Input_matrix_coeff(double A[][N], double, double, double[], int);
void Input_x_p(double, double[], double[], double, int);

//
int main()
{
    setlocale(LC_ALL, "");
    int sz = 9;
    double h = 0.1;
    double x[N] = { NAN }, p[N] = { NAN }, y[N] = { NAN };

    Input_x_p(a, x, p, h, sz);

    double A[N][N];
    zeroing(A);

    Input_matrix_coeff(A, q, h, p, sz);

    double rights[N] = { NAN };

    Input_rights(rights, N, h, p, f, Y0, Yn, sz);

    ShowMatrix(A, x, rights, sz);
    cout << "Значения y на отрезке x = " << x[1] << "..." << x[sz] << " при шаге " << h << endl;
    Thomas(A, rights, y, sz);
    //Proverka(A, y, sz);
}
//
void Thomas(double C[][N], double f[], double y[], int sz)
{
    double S[N] = { NAN };
    double K[N] = { NAN };

    S[1] = -C[0][1] / C[0][0];
    K[1] = f[0] / C[0][0];

    for (int i = 1; i < sz - 1; i++)
    {
        S[i + 1] = -C[i][i + 1] / ((C[i][i - 1] * S[i]) + C[i][i]);
        K[i + 1] = (f[i] - (C[i][i - 1] * K[i])) / ((C[i][i - 1] * S[i]) + C[i][i]);
    }
    y[sz - 1] = ((-C[sz - 1][sz - 2] * K[sz - 1]) - f[sz - 1]) / ((C[sz - 1][sz - 2] * S[sz - 1]) + C[sz - 1][sz - 1]);
    for (int j = sz - 2; j >= 0; j--)
    {
        y[j] = (S[j + 1] * y[j + 1]) + K[j + 1];
    }
    for (int i = 0; i < sz; i++)
    {
        cout << "y[" << "0." << i + 1 << "] = " << y[i] << endl;
    }
}
void zeroing(double Matrix[][N])
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++) Matrix[i][j] = 0;
    }
}
void ShowMatrix(double A[][N], double x[], double rights[], int sz)
{
    cout << "Матрица коэффицентов, A: ";
    int i, j;
    for (i = 0; i < sz; ++i)
    {
        cout << endl;
        for (j = 0; j < sz; ++j) cout << " " << setw(9) << fixed << setprecision(3) << A[i][j];
    }
    //
    cout << endl << endl << "Значения x: " << endl;
    for (int i = 1; i < sz + 1; i++) cout << " " << setw(9) << fixed << setprecision(3) << x[i];
    cout << endl;
    //
    cout << endl << "Правые части уравнений: " << endl;
    for (int i = 0; i < sz; i++) cout << " " << setw(9) << fixed << setprecision(3) << rights[i];
    cout << endl;
}
void Proverka(double matrix[][N], double v[], int sz)
{
    cout << endl << "Проверка: " << endl;
    double vn[N];
    for (int i = 0; i < sz; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < sz; j++) sum = sum + matrix[i][j] * v[j];
        vn[i] = sum;
    }
    for (int i = 0; i < sz; i++)
    {
        cout << vn[i] << endl;
    }
    cout << endl;
}
void Input_rights(double rights[N], int N, double h, double p[], double f, double y0, double yn, int sz)
{
    for (int i = 0; i < sz; i++)
    {
        if (i == 0)
        {
            rights[i] = f - y0 * (1 / (h * h) - p[i + 1] / (2 * h));
        }
        else if (i == sz - 1)
        {
            rights[i] = -(f - yn * (1 / (h * h) + p[i + 1] / (2 * h)));
        }
        else
        {
            rights[i] = f;
        }
    }
}
void Input_matrix_coeff(double A[][N], double q, double h, double p[], int sz)
{
    for (int i = 0; i < sz; i++)
    {
        if (i == 0)
        {
            A[i][i] = (-2 / (h * h) + q);
            A[i][i + 1] = (1 / (h * h) + p[i + 1] / (2 * h));
        }
        else if (i == sz - 1)
        {
            A[i][i - 1] = (1 / (h * h) - p[i + 1] / (2 * h));
            A[i][i] = (-2 / (h * h) + q);
        }
        else
        {
            A[i][i - 1] = (1 / (h * h) - p[i + 1] / (2 * h));
            A[i][i] = (-2 / (h * h) + q);
            A[i][i + 1] = (1 / (h * h) + p[i + 1] / (h * 2));
        }
    }
}
void Input_x_p(double a, double x[N + 1], double p[N + 1], double h, int sz)
{
    for (int i = 0; i < sz + 1; i++)
    {
        x[i] = a;
        p[i] = x[i];
        a = a + h;
    }
}