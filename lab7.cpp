#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace std;

#define N 5
#define MAX_IT 60

double EPS_RES = 1e-8;
double EPS_ERR = 1e-8;

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
            cout << fixed << setprecision(3) << M[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void print_vector(double V[N]) {
    for (int i = 0; i < N; i++) {
        cout << fixed << setprecision(10) << V[i] << " ";
    }
}

void print_info(double v[N], double res, double est, char method) {
    cout << method << " | ";
    print_vector(v);
    cout << "| " << res;
    cout << " | " << est;
    cout << endl;
}

void sum_mm(double M1[N][N], double M2[N][N], double R[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i][j] = M1[i][j] + M2[i][j];
        }
    }
}

void sum_vv(double v1[N], double v2[N], double r[N]) {
    for (int i = 0; i < N; i++) {
        r[i] = v1[i] + v2[i];
    }
}

void sub_vv(double v1[N], double v2[N], double r[N]) {
    for (int i = 0; i < N; i++) {
        r[i] = v1[i] - v2[i];
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

double vector_max_norm(double v[N]) {
    double res = 0.;
    for (int i = 0; i < N; i++) {
        if (abs(v[i]) > abs(res)) {
            res = abs(v[i]);
        }
    }

    return res;
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
    // inicjalizacja zmiennych do metody Jacobiego
    double Dr[N][N] = {};
    double mDr[N][N] = {};
    double LU[N][N] = {};
    for (int i = 0; i < N; i++) {
        Dr[i][i] = 1. / D[i][i];  // D^-1
        mDr[i][i] = -1. / D[i][i];  // -D^-1
    }
    
    sum_mm(L, U, LU);
    mult_mv(Dr, b, Jc);  // C
    mult_mm(mDr, LU, Jm);  // M

    copy(x0, x0 + N, Jx);
    fill(Jxn, Jxn + N, 0.0);
}

bool Jacoby_method(double Jm[][N], double Jc[N], double Jx[N], double Jxn[N]) {
    bool res = false;
    double tmp1[N] = {};
    mult_mv(Jm, Jx, tmp1);
    sum_vv(tmp1, Jc, Jxn);

    // warunki końca iteracji Jacoby
    double diff[N] = {};
    double r[N] = {};
    double tmp2[N] = {};
    sub_vv(Jxn, Jx, diff);
    mult_mv(A, Jxn, tmp2);  
    sub_vv(b, tmp2, r);  

    double norm_diff = vector_max_norm(diff);
    double norm_r = vector_max_norm(r);
    if (norm_r < EPS_RES && norm_diff < EPS_ERR) {
        res = true;
    }

    print_info(Jxn, norm_r, norm_diff, 'J');

    copy(Jxn, Jxn + N, Jx);
    fill(Jxn, Jxn + N, 0.);
    return res;
}

void init_G(double GL[][N], double GU[][N], double Gx[N], double Gxn[N]) {
    // inicjalizacja zmiennych do metody Gaussa-Seidela
    mult_ms(U, -1., GU);  // -U
    sum_mm(L, D, GL);  // L + D

    copy(x0, x0 + N, Gx);
    fill(Gxn, Gxn + N, 0.0);
}

bool Gauss_Seidel_method(double GL[][N], double GU[][N], double Gx[N], double Gxn[N]) {
    bool res = false;
    double tmp[N];
    double p[N];

    // rozwiazanie ukladu rownan
    mult_mv(GU, Gx, tmp);
    sum_vv(tmp, b, p);
    for (int i = 0; i < N; i++) {
        double sum = 0.;
        for (int j = 0; j < i; j++) {
            sum += GL[i][j] * Gxn[j];
        }
        Gxn[i] = (p[i] - sum) / GL[i][i];
    }

    // warunki konca iteracji Gauss-Seidel
    double diff[N] = {};
    double r[N] = {};
    double tmp2[N] = {};
    sub_vv(Gxn, Gx, diff);
    mult_mv(A, Gxn, tmp2);  
    sub_vv(b, tmp2, r);  

    double norm_diff = vector_max_norm(diff);
    double norm_r = vector_max_norm(r);
    if (norm_r < EPS_RES && norm_diff < EPS_ERR) {
        res = true;
    }
    
    print_info(Gxn, norm_r, norm_diff, 'G');

    copy(Gxn, Gxn + N, Gx);
    fill(Gxn, Gxn + N, 0.);
    return res;
}

void init_S(double SL[][N], double SU[N][N], double Sx[N], double Sxn[N]) {
    // inicjalizacja zmiennych do metody SOR
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

bool SOR_method(double SL[][N], double SU[][N], double Sx[N], double Sxn[N]) {
    bool res = false;
    double tmp[N];
    double p[N];

    // rozwiazanie układu rownan
    mult_mv(SU, Sx, tmp);
    sum_vv(tmp, b, p);
    for (int i = 0; i < N; i++) {
        double sum = 0.;
        for (int j = 0; j < i; j++) {
            sum += SL[i][j] * Sxn[j];
        }
        Sxn[i] = (p[i] - sum) / SL[i][i];
    }
    
    // warunki konca iteracji Gauss-Seidel
    double diff[N] = {};
    double r[N] = {};
    double tmp2[N] = {};
    sub_vv(Sxn, Sx, diff);
    mult_mv(A, Sxn, tmp2);  
    sub_vv(b, tmp2, r);  

    double norm_diff = vector_max_norm(diff);
    double norm_r = vector_max_norm(r);
    if (norm_r < EPS_RES && norm_diff < EPS_ERR) {
        res = true;
    }

    print_info(Sxn, norm_r, norm_diff, 'S');

    copy(Sxn, Sxn + N, Sx);
    fill(Sxn, Sxn + N, 0.);
    return res;
}

int main() {
    double Jm[N][N], Jc[N], Jx[N], Jxn[N];  // zmienne metoda Jacobiego
    double GL[N][N], GU[N][N], Gx[N], Gxn[N];  // zmienne metoda Gaussa-Seidela
    double SL[N][N], SU[N][N], Sx[N],  Sxn[N];  // zmienne metoda SOR

    decompose_A();
    init_J(Jm, Jc, Jx, Jxn);
    init_G(GL, GU, Gx, Gxn);
    init_S(SL, SU, Sx, Sxn);

    int i = 0, jr = 0, gr = 0, sr = 0;
    bool jacoby = false;
    bool gauss = false;
    bool sor = false;
    while (true) {
        i++;

        cout << "----------------------------------------- iteracja " << i << " ------------------------------------------" << endl; 
        if (!jacoby) {
            jacoby = Jacoby_method(Jm, Jc, Jx, Jxn);
            jr++;
        }
        if (!gauss) {
            gauss = Gauss_Seidel_method(GL, GU, Gx, Gxn);
            gr++;
        }
        if (!sor) {
            sor = SOR_method(SL, SU, Sx, Sxn);
            sr++;
        }
        cout << endl << endl;

        if (jacoby && gauss && sor) {
            break;
        }

        if (i >= MAX_IT) {
            break;
        }
    }

    cout << "Wyniki końcowe:" << endl;
    cout << "Jacoby       " << " | ";
    print_vector(Jx);
    cout << " | " << jr << " iteracji" << endl; 
    cout << "Gauss-Seidel " << " | ";
    print_vector(Gx);
    cout << " | " << gr << " iteracji" << endl; 
    cout << "SOR          " << " | ";
    print_vector(Sx);
    cout << " | " << sr << " iteracji" << endl; 
}