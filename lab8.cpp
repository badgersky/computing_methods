#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

double x0d = 0.0;
double x1d = M_PI / 4.0;
double x2d = M_PI / 2.0;
long double x0l = 0.0l;
long double x1l = M_PI / 4.0l;
long double x2l = M_PI / 2.0l;

template <typename T>
T forward_2(T x, T h) {
    return (cos(x + h) - cos(x)) / h;
}

template <typename T>
T backward_2(T x, T h) {
    return (cos(x) - cos(x - h)) / h;
}

template <typename T>
T forward_3(T x, T h) {
    return (-3 * cos(x) + 4 * cos(x + h) - cos(x + 2 * h)) / (2 * h);
}

template <typename T>
T backward_3(T x, T h) {
    return (cos(x - 2 * h) - 4 * cos(x - h) + 3 * cos(x)) / (2 * h);
}

template <typename T>
T central_3(T x, T h) {
    return (cos(x + h) - cos(x - h)) / (2 * h);
}


int main() {
    double h1;
    long double h2;

    ofstream file_double("C:\\Users\\szymon\\workbench\\computing_methods\\double_errors.txt");
    ofstream file_long("C:\\Users\\szymon\\workbench\\computing_methods\\longdouble_errors.txt");

    for (h1 = 0.1; h1 > 1e-15; h1 /= 1.1) {
        double err_f2 = abs(forward_2(x0d, h1) - (-sin(x0d)));
        double err_b2 = abs(backward_2(x2d, h1) - (-sin(x2d)));
        double err_f3 = abs(forward_3(x0d, h1) - (-sin(x0d)));
        double err_b3 = abs(backward_3(x2d, h1) - (-sin(x2d)));
        double err_c3 = abs(central_3(x1d, h1) - (-sin(x1d)));

        file_double << h1 << " "
                    << err_f2 << " "
                    << err_b2 << " "
                    << err_f3 << " "
                    << err_b3 << " "
                    << err_c3 << endl;
    }

    for (h2 = 0.1L; h2 > 1e-18L; h2 /= 1.1L) {
        long double err_f2 = fabsl(forward_2(x0l, h2) - (-sinl(x0l)));
        long double err_b2 = fabsl(backward_2(x2l, h2) - (-sinl(x2l)));
        long double err_f3 = fabsl(forward_3(x0l, h2) - (-sinl(x0l)));
        long double err_b3 = fabsl(backward_3(x2l, h2) - (-sinl(x2l)));
        long double err_c3 = fabsl(central_3(x1l, h2) - (-sinl(x1l)));

        file_long << h2 << " "
                  << err_f2 << " "
                  << err_b2 << " "
                  << err_f3 << " "
                  << err_b3 << " "
                  << err_c3 << endl;
    }

    file_double.close();
    file_long.close();

    return 0;
}