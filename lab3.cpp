#include <cmath>
#include <iostream>
#include <iomanip>
#include <random>

using namespace std;

double EPS = 1e-14;
double N = 100;

double func1(double x) {
    return tanh(x) + 2. * (x - 1.);
}

double func2(double x) {
    return sinh(x) + x / 4. - 1.;
}

double dfunc1(double x) {
    return 2. + 1. / (cosh(x)*cosh(x));
}

double dfunc2(double x) {
    return cosh(x) + 1./4.;
}

double bisection_method(double a, double b, double(*func)(double)) {
    double av, bv, x0v, x0;
    int i = 0;
    av = func(a);
    bv = func(b);

    if (av * bv > 0.) {
        throw runtime_error("funkcja nie spełnia założeń!");
    } else {
        while (true) {
            i++;
            x0 = (a + b) / 2.;
            x0v = func(x0);

            cout << endl;
            cout << "---------------- " << i << "-ta iteracja ----------------" << endl;
            cout << "Przybliżenie pierwiastka: " << x0 << endl;
            cout << "Wartość funkcji w x0:     " << x0v << endl;
            cout << "Estymator błędu:          " << (fabs(a - b) / 2.) << endl;

            if ((fabs(a - b) / 2.) < EPS) break;

            if (fabs(x0v) < EPS) break;

            if (i >= N) throw runtime_error("osiągnięto limit iteracji!");

            if (av * x0v < 0.) {b = x0; bv = x0v;}
            else {a = x0; av = x0v;}
        }

        return x0;
    }
}

double secant_method(double a, double b, double(*func)(double)) {
    double x1, x2, x3, av, bv, x1v, x2v, x3v;
    int i = 0;
    av = func(a);
    bv = func(b);
    x1 = a;
    x2 = b;

    if (av * bv > 0.) {
        throw runtime_error("funkcja nie spełnia założeń!");
    } else {
        while (true) {
            i++;
            x1v = func(x1);
            x2v = func(x2);
            
            x3 = x2 - (x2v * (x2 - x1)) / (x2v - x1v);
            x3v = func(x3);

            cout << endl;
            cout << "---------------- " << i << "-ta iteracja ----------------" << endl;
            cout << "Przybliżenie pierwiastka: " << x3 << endl;
            cout << "Wartość funkcji w x3:     " << fabs(x3v) << endl;
            cout << "Estymator błędu:          " << fabs(x2 - x1) << endl;

            if (fabs(x3v) < EPS) break;

            if (fabs(x2 - x3) < EPS) break;

            if (i >= N) throw runtime_error("osiągnięto limit iteracji!");

            x1 = x2;
            x2 = x3;
        }

        return x3;
    }
}

double newton_method(double a, double b, double(*func)(double), double(*dfunc)(double)) {
    double av, bv, xn, xnv, x0, x0v1, x0v2;
    int i = 0;
    av = func(a);
    bv = func(b);
    x0 = (a + b) / 2.;

    if (av * bv > 0.) {
        throw runtime_error("funkcja nie spełnia założeń!");
    } else {
        while (true) {
            i++;
            x0v1 = func(x0);
            x0v2 = dfunc(x0);
            if (x0v2 < EPS) throw runtime_error("pochodna z x0 wynosi 0!");

            xn = x0 - x0v1 / x0v2;
            xnv = func(xn);

            cout << endl;
            cout << "---------------- " << i << "-ta iteracja ----------------" << endl;
            cout << "Przybliżenie pierwiastka: " << xn << endl;
            cout << "Wartość funkcji w x3:     " << fabs(xnv) << endl;
            cout << "Estymator błędu:          " << fabs(x0 - xn) << endl;

            if (fabs(xnv) < EPS) break;

            if (fabs(x0 - xn) < EPS) break;

            if (i >= N) throw runtime_error("osiągnięto limit iteracji!");

            x0 = xn;
        }

        return xn;
    }
}

int main() {
    cout << setprecision(20) << fixed;

    double bires1, bires2, seres1, seres2, neres1, neres2;

    try {
        cout << endl;
        cout << "----------------- METODA BISEKCJI -----------------" << endl;
        bires1 = bisection_method(-10., 10., func1);
        bires2 = bisection_method(-10., 10., func2);
    }
    catch (const runtime_error& e) {
        cerr << "Błąd w metodzie bisekcji: " << e.what() << endl;
    }

    try {
        cout << endl;
        cout << "---------------- METODA SIECZNYCH ----------------" << endl;
        seres1 = secant_method(-10., 10., func1);
        seres2 = secant_method(-10., 10., func2);
    }
    catch (const runtime_error& e) {
        cerr << "Błąd w metodzie siecznych: " << e.what() << endl;
    }

    try {
        cout << endl;
        cout << "----------------- METODA NEWTONA -----------------" << endl;
        neres1 = newton_method(-10., 10., func1, dfunc1);
        neres2 = newton_method(-10., 10., func2, dfunc2);
    }
    catch (const runtime_error& e) {
        cerr << "Błąd w metodzie Newtona: " << e.what() << endl;
    }

    cout << endl;
    cout << "----------------------- \t Porównanie Metod \t -----------------------" << endl;
    cout << "metoda \t\t\t funkcja 1 \t\t\t funkcja 2" << endl;
    cout << "metoda bisekcji\t\t " << bires1 << "\t\t " << bires2 << endl;
    cout << "metoda siecznych\t " << seres1 << "\t\t " << seres2 << endl;
    cout << "metoda Newtona\t\t " << neres1 << "\t\t " << neres2 << endl;
    return 0;
}