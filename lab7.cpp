#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace std;

#define N 5

double A[N][N] = {
    {50, 5, 4, 3, 2},
    {1, 40, 1, 2, 3},
    {4, 5, 30, -5, -4},
    {-3, -2, -1, 20, 0},
    {1, 2, 3, 4, 30}
};

double b[N] = {140, 67, 62, 89, 153};

double x0[N]={6, 6, 6, 6, 6};

double L[N][N] = {}, D[N][N] = {}, U[N][N] = {};

double o = 0.5;

void print_matrix(double M[][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << setw(10) << fixed << setprecision(2) << M[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void print_vector(double V[N]) {
    for (int i = 0; i < N; i++) {
        cout << setw(10) << fixed << setprecision(2) << V[i] << " ";
    }
    cout << endl;
}

void sum_mm(double M1[N][N], double M2[N][N], double R[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i][j] = M1[i][j] + M2[i][j];
        }
    }
}

void mult_mm(double B1[][N], double B2[][N], double R[][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i][j] = 0.0;
            for (int k = 0; k < N; k++) {
                R[i][j] += B1[i][k] * B2[k][j];
            }
        }
    }
}

void mult_ms(double B1[N][N], double s, double R[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i][j] = s * B1[i][j];
        }
    }
}

void mult_mv(double B[][N], double v[N], double r[N]) {
    for (int i = 0; i < N; i++) {
        r[i] = 0.0;
        for (int j = 0; j < N; j++) {
            r[i] += B[i][j] * v[j];
        }
    }
}

void decompose_A() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) D[i][j] = A[i][j];
            else if (i > j) L[i][j] = A[i][j];
            else U[i][j] = A[i][j];
        }
    }
}

void init_J(double Jm[][N], double Jc[N], double Jx[N], double Jxn[N]) {
    double Dr[N][N] = {};
    double mDr[N][N] = {};
    double LU[N][N] = {};
    for (int i = 0; i < N; i++) {
        Dr[i][i] = 1. / D[i][i];
        mDr[i][i] = -1. / D[i][i];
    }
    
    sum_mm(L, U, LU);
    mult_mv(Dr, b, Jc);  // wektor C
    mult_mm(mDr, LU, Jm);  // macierz M

    copy(x0, x0 + N, Jx);
    fill(Jxn, Jxn + N, 0.0);
}

void init_G(double GL[][N], double GU[][N], double Gx[N], double Gxn[N]) {
    mult_ms(U, -1., GU);
    sum_mm(L, D, GL);  // L + D

    copy(x0, x0 + N, Gx);
    fill(Gxn, Gxn + N, 0.0);
}

void init_S(double SL[][N], double SU[N][N], double Sx[N], double Sxn[N]) {
    // L + 1/o * D
    double oD1[N][N] = {};
    mult_ms(D, (1. / o), oD1);
    sum_mm(L, oD1, SL);

    // - [(1 - o) * D + U]
    double Su[N][N] = {};
    double oD2[N][N] = {};
    mult_ms(D, (1. - (1. / o)), oD2);
    sum_mm(oD2, U, Su);
    mult_ms(Su, -1, SU);

    copy(x0, x0 + N, Sx);
    fill(Sxn, Sxn + N, 0.0);
}

void Jacoby_method(double Jm[][N], double Jc[N], double Jx[N], double Jxn[N]) {
    cout << "Jacoby" << endl;
}

void Gauss_Seidel_method(double Gm[][N], double Gc[N], double Gx[N], double Gxn[N]) {
    cout << "Gauss-Seidel" << endl;
}

void SOR_method(double Sm[][N], double Sc[N], double Sx[N], double Sxn[N]) {
    cout << "SOR" << endl;
}

int main() {
    double Jm[N][N], GL[N][N], GU[N][N], SL[N][N], SU[N][N];
    double Jc[N], Jx[N], Gx[N], Sx[N], Jxn[N], Gxn[N], Sxn[N];

    decompose_A();
}