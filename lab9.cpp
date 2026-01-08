#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

#define N 1000
#define M 69
#define EPS1 1e-13
#define EPS2 1e-13

using namespace std;

void init(double a, double h, double A, double B, double C, double lower[], double upper[], double diag[], double r[]) {
    double x;

    for (int i = 1; i < N - 1; i++) {
        lower[i - 1] = A;
        diag[i] = B;
        upper[i] = C;

        x = a + i * h;
        r[i] = - x * x * x / 2.;
    }

    diag[0] = 1.;
    r[0] = 2.;
    upper[0] = 0.;

    diag[N - 1] = 1.;
    r[N - 1] = -2.;
    lower[N - 2] = 0.;
}

void thomas_algorithm(double lower[], double diag[], double upper[], double r[], double x[]) {
    double c_star[N - 1];
    double d_star[N];

    c_star[0] = upper[0] / diag[0];
    d_star[0] = r[0] / diag[0];

    for (int i = 1; i < N; i++) {
        double m = diag[i] - lower[i - 1] * c_star[i - 1];
        if (i < N - 1) {
            c_star[i] = upper[i] / m;
        }
        d_star[i] = (r[i] - lower[i - 1] * d_star[i - 1]) / m;
    }

    x[N - 1] = d_star[N - 1];
    for (int i = N - 2; i >= 0; i--) {
        x[i] = d_star[i] - c_star[i] * x[i + 1];
    }
}

double get_exact_u(double x) {
    double s5 = sqrt(5.0);
    double exp1 = exp((-1 - s5) * (-1 + x));
    double exp2 = exp((-1 + s5) * x);
    double exp3 = exp(1 + s5 + (-1 + s5) * x);
    double exp4 = exp(2 * s5 - (1 + s5) * x);
    double exp5 = exp(2 * s5);

    double r1 = 9 - 95 * exp1 + 55 * exp2 + 95 * exp3 - 55 * exp4 + 2 * x * (6 + x * (3 + 2 * x)) - exp5 * (9 + 2 * x * (6 + x * (3 + 2 * x)));
    double r2 = -1 + 1.0/tanh(s5);

    return - r1 * r2 / 64.0;
}

double bisection_method(double a, double b, double(*func)(double)) {
    double av, bv, xnv, xn;
    int i = 0;
    av = func(a);
    bv = func(b);

    if (av * bv > 0.) {
        throw runtime_error("funkcja nie spełnia założeń!");
    } else {
        while (true) {
            i++;
            xn = (a + b) / 2.;
            xnv = func(xn);

            cout << endl;
            cout << "---------------- " << i << "-ta iteracja ----------------" << endl;
            cout << "Przybliżenie pierwiastka xn: " << xn << endl;
            cout << "Wartość funkcji w xn:        " << fabs(xnv) << endl;
            cout << "Estymator błędu:             " << (fabs(a - b) / 2.) << endl;

            if ((fabs(a - b) / 2.) < EPS1 && fabs(xnv) < EPS2) break;

            if (i >= N) throw runtime_error("osiągnięto limit iteracji!");

            if (av * xnv < 0.) {b = xn; bv = xnv;}
            else {a = xn; av = xnv;}
        }

        return xn;
    }
}

void shooting_method(double a, double h, double u_sh[], double p, double A, double B, double C) {
    double x, f;
    
    u_sh[0] = 2;
    u_sh[1] = u_sh[0] + h * p;

    for (int i = 2; i < N; i++) {
        x = a + i * h;
        f = - x * x * x / 2.;
        u_sh[i] = (-1. * A * u_sh[i - 2] + -1. * B * u_sh[i - 1] - f) / C;
    }
}

void get_exact_solution(double a, double h, double u_exact[]) {
    double x;
    for (int i = 0; i < N; i++) {
        x = a + i * h;
        u_exact[i] = get_exact_u(x);
    }
}

int main() {
    double a = 0.;
    double b = 1.;
    double h = abs(a - b) / (N - 1);
    double A = 1. / (h * h) - 1. / h;
    double B = -2. / (h * h) - 4.;
    double C = 1. / (h * h) + 1. / h;

    double lower[N - 1];
    double upper[N - 1];
    double diag[N];
    double r[N];
    double u[N];
    double u_exact[N];
    double u_sh[N];

    init(a, h, A, B, C, lower, upper, diag, r);
    thomas_algorithm(lower, diag, upper, r, u);
    get_exact_solution(a, h, u_exact);
    shooting_method(a, h, u_sh, 0., A, B, C);

    for (int i = 0; i < N; i++) {
        cout << setw(10) << setprecision(12) << fixed << u_sh[i] << "   " << u_exact[i] << "   " << abs(u_exact[i] - u_sh[i]) << endl;
    }
}