// ��������� ������� ������� ����������������� ��������� ���� y'' + p(x)y' + q(x)y = f(x) � ������ y'' + xy' - 2y = -1.5
#include <iostream>
#include <stdlib.h>
#include <iomanip>
const int N = 9; 
double h = 0.1;  
double a = 0; 
double b = 1; 
double q = -2; 
double f = -1.5; 
#include "Class_Matrix_Vector.h"
using namespace std;
void init_matrix_coeff(Matrix&, Vector&);
void init_x_p(double a, Vector&, Vector&, double);
void init_rights(Vector& rights, double h, Vector& p, double f, const double y0, const double yn);
void Thomas_for_classes(Matrix& C, Vector& f, Vector& y);

int main()
{
    setlocale(LC_ALL, "");
    const double y0 = 2;
    const double yn = 2 * sqrt(2); 
    Vector vx(N), vp(N), vy(N);
    init_x_p(a, vx, vp, h);
    Matrix A(N, N);
    init_matrix_coeff(A, vp);
    cout << "������� ������������";
    A.Show();
    Vector v_rights(N);
    init_rights(v_rights, h, vp, f, y0, yn);
    cout << endl << endl << "������ ����� ���������";
    v_rights.Show();
    cout << endl << "�������� y �� ������� x = " << vx(1) << "..." << vx(N) << " ��� ���� " << h << endl;
    Thomas_for_classes(A, v_rights, vy);
}
 
    
    
    