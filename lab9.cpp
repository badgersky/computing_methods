#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

#define N 1000
#define M 100
#define EPS1 1e-13
#define EPS2 1e-13

using namespace std;

void init(double a, double h, double A, double B, double C, double lower[], double upper[], double diag[], double r[], int n) {
    double x;

    for (int i = 1; i < n - 1; i++) {
        lower[i - 1] = A;
        diag[i] = B;
        upper[i] = C;

        x = a + i * h;
        r[i] = - x * x * x / 2.;
    }

    diag[0] = 1.;
    r[0] = 2.;
    upper[0] = 0.;

    diag[n - 1] = 1.;
    r[n - 1] = -2.;
    lower[n - 2] = 0.;
}

void thomas_algorithm(double lower[], double diag[], double upper[], double r[], double x[], int n) {
    double c_star[n - 1];
    double d_star[n];

    c_star[0] = upper[0] / diag[0];
    d_star[0] = r[0] / diag[0];

    for (int i = 1; i < n; i++) {
        double m = diag[i] - lower[i - 1] * c_star[i - 1];
        if (i < n - 1) {
            c_star[i] = upper[i] / m;
        }
        d_star[i] = (r[i] - lower[i - 1] * d_star[i - 1]) / m;
    }

    x[n - 1] = d_star[n - 1];
    for (int i = n - 2; i >= 0; i--) {
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

void get_exact_solution(double a, double h, double u_exact[], int n) {
    double x;
    for (int i = 0; i < n; i++) {
        x = a + i * h;
        u_exact[i] = get_exact_u(x);
    }
}

void shoot(double p, double h, double u_sh[], double A, double B, double C, int n) {
    u_sh[0] = 2.0;
    u_sh[1] = u_sh[0] + h * p;

    for(int i = 2; i < n; i++) {
        double x = i * h;
        double f = - x*x*x / 2.0;
        u_sh[i] = (f - A * u_sh[i - 2] - B * u_sh[i - 1]) / C;
    }
}

double shooting_function(double p, double h, double u_sh[], double A, double B, double C, int n) {
    shoot(p, h, u_sh, A, B, C, n);
    return u_sh[n - 1] - (-2.0);
}

void shooting_method(double h, double u_sh[], double A, double B, double C, int n) {
    double p1 = -10.0;
    double p2 = 10.0;
    double u_temp[n];

    double f1 = shooting_function(p1, h, u_temp, A, B, C, n);
    double f2 = shooting_function(p2, h, u_temp, A, B, C, n);

    if(f1 * f2 > 0) {
        exit(1);
    }

    double s_star;

    for(int i = 0; i < M; i++) {
        s_star = 0.5 * (p1 + p2);
        double f_mid = shooting_function(s_star, h, u_sh, A, B, C, n);

        if(f_mid * f1 < 0) {
            p2 = s_star;
            f2 = f_mid;
        } else {
            p1 = s_star;
            f1 = f_mid;
        }

        if(fabs(p2 - p1) / 2. < EPS1 && fabs(f_mid) < EPS2) break;
    }
}

int main() {
    for (int n = 10; n <= 1000; n += 10) {
        double a = 0.;
        double b = 1.;
        double h = abs(a - b) / (n - 1);
        double A = 1. / (h * h) - 1. / h;
        double B = -2. / (h * h) - 4.;
        double C = 1. / (h * h) + 1. / h;
        double p_final;
        double p1 = -10., p2 = -9.;

        double lower[n - 1];
        double upper[n - 1];
        double diag[n];
        double r[n];
        double u_thomas[n];
        double u_exact[n];
        double u_shooting[n];

        init(a, h, A, B, C, lower, upper, diag, r, n);
        thomas_algorithm(lower, diag, upper, r, u_thomas, n);
        get_exact_solution(a, h, u_exact, n);
        shooting_method(h, u_shooting, A, B, C, n);
    }
}