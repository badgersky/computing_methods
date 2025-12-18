#include <iostream>
#include <algorithm>
#include <iomanip>
#include <cmath>

using namespace std;

#define N 5 

double A[N][N] = {
    {5., 4., 3., 2., 1.},
    {10., 8., 7., 6., 5.},
    {-1., 2., -3., 4., -5.},
    {6., 5., -4., 3., -2.},
    {1., 2., 3., 4., 5.}
};

double b[N] = {37., 99., -9., 12., 53.};

int indexes[N] = {0, 1, 2, 3, 4};

double eps = 1e-14;

void print_matrix(double A[][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << setw(10) << fixed << setprecision(5) << A[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void print_L(double A[][N]) {
    cout << "Macierz L:" << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i > j)
                cout << setw(10) << fixed << setprecision(5) << A[i][j] << " ";
            else if (i == j)
                cout << setw(10) << fixed << setprecision(5) << 1.0 << " ";
            else
                cout << setw(10) << fixed << setprecision(5) << 0.0 << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void print_U(double A[][N]) {
    cout << "Macierz U:" << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i <= j)
                cout << setw(10) << fixed << setprecision(5) << A[i][j] << " ";
            else
                cout << setw(10) << fixed << setprecision(5) << 0.0 << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void print_vector(double x[], char c) {
    cout << "Rozwiązanie dla wektora " << c << ":" << endl;
    for (int i = 0; i < N; i++)
        cout << fixed << setprecision(5) << c << i + 1 << " = " << x[i] << endl;
    cout << endl;
}

void decomposition(double A[][N], int indexes[N]) {
    double p, w;
    cout << "------------------ macierz pierwotna -------------------" << endl;
    print_matrix(A);
    for (int i = 0; i < N; i++) {
        p = A[i][i];
        
        cout << "---------------------- " << "iteracja " << i + 1 << " ----------------------" << endl;
        if (abs(p) < eps) {
            int max_row = i;
            for (int j = i + 1; j < N; j++) {
                if (abs(A[j][i]) > abs(A[max_row][i])) {
                    max_row = j;
                }
            }

            if (max_row != i) {
                for (int j = 0; j < N; j++) {
                    swap(A[i][j], A[max_row][j]);
                }
                swap(indexes[i], indexes[max_row]);
                p = A[i][i];
            }

            if (abs(p) < eps) {
                return;
            }
        }
        
        cout << "Element podstawowy: " << p << endl;
        for (int j = i + 1; j < N; j++) {
            w = A[j][i] / p;
            for (int k = i + 1; k < N; k++) {
                A[j][k] -= w * A[i][k];
            }
            A[j][i] = w;
        }
        print_matrix(A);
    }
}

void solve(double A[][N], double b[], int indexes[]) {
    double y[N];
    double x[N];

    for (int i = 0; i < N; i++) {
        y[i] = b[indexes[i]];
        for (int j = 0; j < i; j++)
            y[i] -= A[i][j] * y[j];
    }
    print_vector(y, 'y');

    for (int i = N - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < N; j++)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
    print_vector(x, 'x');
}

int main() {
    decomposition(A, indexes);
    print_U(A);
    print_L(A);
    solve(A, b, indexes);
}
