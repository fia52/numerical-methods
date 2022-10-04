#pragma once
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
using namespace std;
class Array
{
public:
    Array();
    Array(int);
    Array(int, int);
    ~Array();
    void grow(int);
    virtual void Show() = 0;
protected:
    void init(const double*, int, int);
    double* p;
    int size;
    int error;
};
//  protected function init
void Array::init(const double* pa, int m, int n)
{
    int i, j;
    size = m * n;
    p = new double[size];
    assert(p != 0);
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < n; ++j)
            p[i * n + j] = (pa != 0) ? pa[i * n + j] : 0.0;
    }
}
// public function grow
void Array::grow(int factor)
{
    double* oldp = p;
    int  oldsize = size;
    size = size * factor;
    p = new double[size];
    for (int i = 0; i < oldsize; ++i)
        p[i] = oldp[i];
    for (int i = oldsize; i < size; ++i)
        p[i] = 0.0;
    delete[] oldp;
}
// implicit constructor
Array::Array(void)
{
    init(0, 0, 0);
}
// 1 constructor
Array::Array(int m)
{
    if (m < 0)
    {
        cout << "\nNegative size of vector !";
        exit(EXIT_FAILURE);
    }
    init(0, m, 1);
}
// 2 constructor
Array::Array(int m, int n)
{
    if (m < 0 || n < 0)
    {
        cerr << "\nNegative size of matrix !";
        exit(EXIT_FAILURE);
    }
    init(0, m, n);
}
// destructor
Array::~Array()
{
    delete[] p;
}
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///                                                                           ///
///                         DERIVED CLASS VECTOR                              ///
///                                                                           ///
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
class Vector :public Array
{
public:
    Vector() :Array() {}
    Vector(int m) :Array(m)
    {
        row = m;
    }
    Vector(const Vector& v)
    {
        init(v.p, v.row, 1);
    }
    Vector(double* v, int sz)
    {
        init(v, sz, 1);
    }
    int get_size() const
    {
        return row;
    }
    friend double Norm(Vector& v)
    {
        int i;
        double nrm = 0.0;
        for (i = 1; i <= v.row; ++i) nrm = nrm + v(i) * v(i);
        nrm = sqrt(nrm);
        return nrm;
    }
    friend double MAXVAL(const Vector& v)
    {
        int i;
        double max = v(1);
        for (i = 1; i <= v.row; ++i) max = fmax(max, v(i));
        return max;
    }
    double& operator()(int i) const
    {
        assert(i >= 1 && i <= size);
        return p[i - 1];
    }
    Vector& operator=(const Vector& v)
    {
        if (this == &v)return *this;
        delete p;
        init(v.p, v.row, 1);
        return *this;
    }
    Vector& operator=(const double var)
    {
        int i;
        for (i = 0; i <= row; ++i) p[i] = var;
        return *this;
    }
    friend Vector operator-(const Vector& v1, const Vector& v2)
    {
        int i, ia, ib;
        ia = v1.row;
        ib = v2.row;
        Vector vc(ia);
        if (ia == ib)
        {
            for (i = 1; i <= ia; ++i)
                vc(i) = v1(i) - v2(i);
        }
        else
        {
            cerr << "\nBad in vector addition !";
            exit(EXIT_FAILURE);
        }
        return vc;
    }
    friend Vector operator+(const Vector& v1, const Vector& v2)
    {
        int i, ia, ib;
        ia = v1.row;
        ib = v2.row;
        Vector vc(ia);
        if (ia == ib)
        {
            for (i = 1; i <= ia; ++i)
                vc(i) = v1(i) + v2(i);
        }
        else
        {
            cerr << "\nBad in vector addition !";
            exit(EXIT_FAILURE);
        }
        return vc;
    }
    friend Vector operator+=(Vector& v1, const Vector& v2)
    {
        int i, ia, ib;
        ia = v1.row;
        ib = v2.row;
        if (ia == ib)
        {
            for (i = 1; i <= ia; ++i)
                v1(i) = v1(i) + v2(i);
        }
        else
        {
            cerr << "\nBad in vector addition!";
            exit(EXIT_FAILURE);
        }
        return v1;
    }
    friend Vector operator-=(Vector& v1, const Vector& v2)
    {
        int i, ia, ib;
        ia = v1.row;
        ib = v2.row;
        if (ia == ib)
        {
            for (i = 1; i <= ia; ++i)
                v1(i) = v1(i) - v2(i);
        }
        else
        {
            cerr << "\nBad in vectors addition!";
            exit(EXIT_FAILURE);
        }
        return v1;
    }
    friend Vector& operator*(Vector& v, double var)
    {
        int i;
        if (var <= DBL_EPSILON)
        {
            for (i = 1; i <= v.row; ++i)
                v(i) = 0.0;
        }
        else
        {
            for (i = 1; i <= v.row; ++i)
                v(i) = v(i) * var;
        }
        return v;
    }
    friend Vector& operator*(double var, Vector& v)
    {
        int i;
        if (var <= DBL_EPSILON)
        {
            for (i = 1; i <= v.row; ++i)
                v(i) = 0.0;
        }
        else
        {
            for (i = 1; i <= v.row; ++i)
                v(i) = v(i) * var;
        }
        return v;
    }
    void Show()
    {
        cout << "\nVector values :" << endl;
        for (int i = 0; i < size; ++i)
            cout << p[i] << endl;
    }
private:
    int row;
};
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///                                                                           ///
///                         DERIVED CLASS MATRIX                              ///
///                                                                           ///
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
class Matrix :public Array
{
public:
    Matrix() :Array() {}
    Matrix(int m, int n) :Array(m, n)
    {
        row = m;
        col = n;
    }
    Matrix(const Matrix& m)
    {
        init(m.p, m.row, m.col);
    }
    Matrix(double m[N][N], int sz1, int sz2) : Matrix(sz1, sz2)
    {
        int i, j;
        size = sz1 * sz2;
        p = new double[size];
        for (i = 0; i < sz1; ++i)
        {
            for (j = 0; j < sz2; ++j)
                p[i * sz1 + j] = (m[i][j] != 0) ? m[i][j] : 0.0;
        }
    }
    Matrix& operator=(const Matrix& m)
    {
        if (this == &m) return *this;
        delete p;
        init(m.p, m.row, m.col);
        return *this;
    }
    Matrix& operator=(const double var)
    {
        int i, j;
        for (i = 0; i < row; ++i) for (j = 0; j < col; ++j) p[i * col + j] = var;
        return *this;
    }
    double& operator() (int i, int j) const
    {
        assert((i >= 1 && i <= row) && (j >= 1 && j <= col));
        return p[(i - 1) * col + (j - 1)];
    }
    Matrix Transpose()
    {
        int k;
        double* oldp = p;
        p = new double[size];
        k = row; row = col; col = k;
        for (int i = 0; i < row; ++i)for (int j = 0; j < col; ++j)
            p[i * col + j] = oldp[j * row + i];
        delete[] oldp;
        return *this;
    }
    friend Matrix ABS(const Matrix& m)
    {
        Matrix mabs(m.row, m.col);
        int i, j;
        for (i = 1; i <= m.row; ++i) for (j = 1; j <= m.col; ++j) mabs(i, j) = abs(m(i, j));
        return mabs;
    }
    int get_rows() const
    {
        return row;
    }
    int get_cols() const
    {
        return col;
    }
    int get_size() const
    {
        return row * col;
    }
    friend Matrix operator+(const Matrix& m1, const Matrix& m2)
    {
        int i, j, ia, ja, ib, jb;
        ia = m1.row;
        ja = m1.col;
        ib = m2.row;
        jb = m2.col;
        Matrix mn(ib, jb);
        if ((ia == ib) && (ja == jb))
        {
            for (i = 1; i <= jb; ++i)
            {
                for (j = 1; j <= jb; ++j)
                    mn(i, j) = m1(i, j) + m2(i, j);
            }
        }
        else
        {
            cerr << "\nIncompability in matrix addition!";
            exit(EXIT_FAILURE);
        }
        return mn;
    }
    friend Matrix operator-(const Matrix& m1, const Matrix& m2)
    {
        int i, j, ia, ja, ib, jb;
        ia = m1.row;
        ja = m1.col;
        ib = m2.row;
        jb = m2.col;
        Matrix mn(ib, jb);
        if ((ia == ib) && (ja == jb))
        {
            for (i = 1; i <= jb; ++i)
            {
                for (j = 1; j <= jb; ++j)
                    mn(i, j) = m1(i, j) - m2(i, j);
            }
        }
        else
        {
            cerr << "\nIncompability in matrix addition!";
            exit(EXIT_FAILURE);
        }
        return mn;
    }
    friend Vector operator*(const Matrix& m, const Vector& v)
    {
        int sz = m.row;
        int i, j, jm, iv;
        double sum;
        jm = m.col; iv = v.get_size();
        if (jm != iv)
        {
            cerr << "\nIncompability in matrix to vector multiplication!";
            exit(EXIT_FAILURE);
        }
        //
        Vector vn(sz);
        for (i = 1; i <= m.row; ++i)
        {
            sum = 0.0;
            for (j = 1; j <= m.col; ++j) sum = sum + m(i, j) * v(j);
            vn(i) = sum;
        }
        return vn;
    }
    friend Matrix operator*(const Matrix& m1, const Matrix& m2)
    {
        int i, j, k, ia, ib, ja, jb;
        double sum;
        ia = m1.row;
        ja = m1.col;
        ib = m2.row;
        jb = m2.col;
        Matrix mn(ia, jb);
        if (ja == ib)
        {
            for (i = 1; i <= ia; ++i)
            {
                for (j = 1; j <= ib; ++j)
                {
                    sum = 0.0;
                    for (k = 1; k <= ja; ++k)
                    {
                        sum += m1(j, k) * m2(k, j);
                    }
                    mn(i, j) = sum;
                }
            }
        }
        else
        {
            cerr << "\nIncompability in matrix multiplication!";
            exit(EXIT_FAILURE);
        }
        return mn;
    }
    friend Matrix operator*=(Matrix& m1, const Matrix& m2)
    {
        int i, j, k, ia, ib, ja, jb;
        double sum;
        ia = m1.row;
        ja = m1.col;
        ib = m2.row;
        jb = m2.col;
        Matrix mn(ia, ja);
        if (ja == ib)
        {
            for (i = 1; i <= ia; ++i)
            {
                for (j = 1; j <= ib; ++j)
                {
                    sum = 0.0;
                    for (k = 1; k <= ja; ++k)
                    {
                        sum += m1(j, k) * m2(k, j);
                    }
                    mn(i, j) = sum;
                }
            }
        }
        else
        {
            cerr << "\nIncompability in matrix multiplication!";
            exit(EXIT_FAILURE);
        }
        m1 = mn;
        //mn.grow(0);
        return m1;
    }
    friend Matrix operator*(const Matrix& m, double var)
    {
        int i, j;
        if (var <= DBL_EPSILON)
        {
            for (i = 1; i <= m.row; ++i)
            {
                for (j = 1; j <= m.col; ++j)
                    m(i, j) = 0.0;
            }
        }
        else
        {
            for (i = 1; i <= m.row; ++i)
            {
                for (j = 1; j <= m.col; ++j)
                    m(i, j) = m(i, j) * var;
            }
        }
        return m;
    }
    friend Matrix operator*(double var, const Matrix& m)
    {
        int i, j;
        if (var <= DBL_EPSILON)
        {
            for (i = 1; i <= m.row; ++i)
            {
                for (j = 1; j <= m.col; ++j)
                    m(i, j) = 0.0;
            }
        }
        else
        {
            for (i = 1; i <= m.row; ++i)
            {
                for (j = 1; j <= m.col; ++j)
                    m(i, j) = var * m(i, j);
            }
        }
        return m;
    }
    friend Matrix operator+=(Matrix& m1, const Matrix& m2)
    {
        int i, j, ia, ja, ib, jb;
        ia = m1.row;
        ja = m1.col;
        ib = m2.row;
        jb = m2.col;
        if ((ia == ib) && (ja == jb))
        {
            for (i = 1; i <= jb; ++i)
            {
                for (j = 1; j <= jb; ++j)
                    m1(i, j) = m1(i, j) + m2(i, j);
            }
        }
        else
        {
            cerr << "\nIncompability in matrix addition!";
            exit(EXIT_FAILURE);
        }
        return m1;
    }
    friend Matrix operator-=(Matrix& m1, const Matrix& m2)
    {
        int i, j, ia, ja, ib, jb;
        ia = m1.row;
        ja = m1.col;
        ib = m2.row;
        jb = m2.col;
        if ((ia == ib) && (ja == jb))
        {
            for (i = 1; i <= jb; ++i)
            {
                for (j = 1; j <= jb; ++j)
                    m1(i, j) = m1(i, j) - m2(i, j);
            }
        }
        else
        {
            cerr << "\nIncompability in matrix addition!";
            exit(EXIT_FAILURE);
        }
        return m1;
    }
    void Show()
    {
        int i, j;
        cout << "\nMatrix elements:";
        for (i = 0; i < row; ++i)
        {
            cout << endl;
            for (j = 0; j < col; ++j)
                cout << " " << setw(9) << fixed << setprecision(3) << p[i * col + j];
        }
    }
private:
    int row;
    int col;
};
void addtwo(Vector&, int);
void addtwo(Vector& v, int size)
{
    for (int i = 1; i <= size; ++i)
        v(i) = v(i) + 2.;
}

void init_matrix_coeff(Matrix& M, Vector& p)
{
    for (int i = 1; i <= N; i++)
    {
        if (i == 1)
        {
            M(i, i) = (-2 / (h * h) + q);
            M(i, i + 1) = (1 / (h * h) + p(i) / (2 * h));
        }
        else if (i == N)
        {
            M(i, i - 1) = (1 / (h * h) - p(i) / (2 * h));
            M(i, i) = (-2 / (h * h) + q);
        }
        else
        {
            M(i, i - 1) = (1 / (h * h) - p(i) / (2 * h));
            M(i, i) = (-2 / (h * h) + q);
            M(i, i + 1) = (1 / (h * h) + p(i) / (h * 2));
        }
    }
}
void init_x_p(double a, Vector& x, Vector& p, double h)
{
    for (int i = 1; i <= N; i++)
    {
        x(i) = a + h;
        p(i) = x(i);
        a = a + h;
    }
}
void init_rights(Vector& rights, double h, Vector& p, double f, const double y0, const double yn)
{
    for (int i = 1; i <= N; i++)
    {
        if (i == 1)
        {
            rights(i) = f - y0 * (1 / (h * h) - p(i) / (2 * h));
        }
        else if (i == N)
        {
            rights(i) = -(f - yn * (1 / (h * h) + p(i) / (2 * h)));
        }
        else
        {
            rights(i) = f;
        }
    }
}
void Thomas_for_classes(Matrix& C, Vector& f, Vector& y)
{
    Vector S(N), K(N);
   
    S(2) = -C(1,2) / C(1,1);
    K(2) = f(1) / C(1,1);

    for (int i = 2; i < N; i++)
    {
        S(i + 1) = -C(i, i + 1) / ((C(i, i - 1) * S(i)) + C(i, i));
        K(i + 1) = (f(i) - (C(i, i - 1) * K(i))) / ((C(i, i - 1) * S(i)) + C(i, i));
    }
    y(N) = ((-C(N, N - 1) * K(N)) - f(N)) / ((C(N, N - 1) * S(N)) + C(N, N));
    for (int j = N - 1; j >= 1; j--)
    {
        y(j) = (S(j + 1) * y(j + 1)) + K(j + 1);
    }
    for (int i = 1; i <= N; i++)
    {
        cout << "y[" << "0." << i << "] = " << y(i) << endl;
    }
}