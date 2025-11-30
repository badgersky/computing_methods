#include <cmath>
#include <iostream>
#include <iomanip>
#include <random>

using namespace std;

double EPS1 = 1e-14;
double EPS2 = 1e-14;
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

double phi_func1(double x) {
    return 1. - tanh(x) / 2.;
}

double phi_func2(double x) {
    return asinh(1. - x / 4.);
}

double phi_func1_d(double x) {
    return (-1. / 2) * (1. / (cosh(x) * cosh(x)));
}

double phi_func2_d(double x) {
    return 1. / sqrt(x * x - 8. * x + 32.);
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

double secant_method(double a, double b, double(*func)(double)) {
    double x1, x2, xn, av, bv, x1v, x2v, xnv;
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
            
            xn = x2 - (x2v * (x2 - x1)) / (x2v - x1v);
            xnv = func(xn);

            cout << endl;
            cout << "---------------- " << i << "-ta iteracja ----------------" << endl;
            cout << "Przybliżenie pierwiastka xn: " << xn << endl;
            cout << "Wartość funkcji w xn:        " << fabs(xnv) << endl;
            cout << "Estymator błędu:             " << fabs(x2 - xn) << endl;

            if (fabs(xnv) < EPS1 && fabs(x2 - xn) < EPS2) break;

            if (i >= N) throw runtime_error("osiągnięto limit iteracji!");

            x1 = x2;
            x2 = xn;
        }

        return xn;
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
            if (x0v2 < EPS1) throw runtime_error("pochodna z x0 wynosi 0!");

            xn = x0 - x0v1 / x0v2;
            xnv = func(xn);

            cout << endl;
            cout << "---------------- " << i << "-ta iteracja ----------------" << endl;
            cout << "Przybliżenie pierwiastka xn: " << xn << endl;
            cout << "Wartość funkcji w xn:        " << fabs(xnv) << endl;
            cout << "Estymator błędu:             " << fabs(x0 - xn) << endl;

            if (fabs(xnv) < EPS1 && fabs(x0 - xn) < EPS2) break;

            if (i >= N) throw runtime_error("osiągnięto limit iteracji!");

            x0 = xn;
        }

        return xn;
    }
}

double picard_method(double a, double b, double(*func)(double), double(*phi_func)(double), double(*dphi_func)(double)) {
    double x0, xn, av, bv, xnv;
    int i = 0;
    av = func(a);
    bv = func(b);
    x0 = (a + b) / 2.;

    if (av * bv > 0.) {
        throw runtime_error("funkcja nie spełnia założeń!");
    } else if (fabs(dphi_func(x0)) < 0) {
        throw runtime_error("funkcja nie jest zbieżna");
    } else {
        while (true) {
            i++;
            xn = phi_func(x0);
            xnv = func(xn);
            
            cout << endl;
            cout << "---------------- " << i << "-ta iteracja ----------------" << endl;
            cout << "Przybliżenie pierwiastka xn: " << xn << endl;
            cout << "Wartość funkcji w xn:        " << fabs(xnv) << endl;
            cout << "Estymator błędu:             " << fabs(x0 - xn) << endl;

            if (fabs(xnv) < EPS1 && fabs(x0 - xn) < EPS2) break;

            if (i >= N) throw runtime_error("osiągnięto limit iteracji!");

            x0 = xn;
        }

        return xn;
    }
}

int main() {
    cout << setprecision(20) << fixed;

    double bires1, bires2, seres1, seres2, neres1, neres2, pires1, pires2;

    try {
        cout << endl;
        cout << "----------------- METODA BISEKCJI 1 -----------------" << endl;
        bires1 = bisection_method(-10., 10., func1);
        cout << endl;
        cout << "----------------- METODA BISEKCJI 2 -----------------" << endl;
        bires2 = bisection_method(-10., 10., func2);
    }
    catch (const runtime_error& e) {
        cerr << "Błąd w metodzie bisekcji: " << e.what() << endl;
    }

    try {
        cout << endl;
        cout << "---------------- METODA SIECZNYCH 1 ----------------" << endl;
        seres1 = secant_method(-10., 10., func1);
        cout << endl;
        cout << "---------------- METODA SIECZNYCH 2 ----------------" << endl;
        seres2 = secant_method(-10., 10., func2);
    }
    catch (const runtime_error& e) {
        cerr << "Błąd w metodzie siecznych: " << e.what() << endl;
    }

    try {
        cout << endl;
        cout << "----------------- METODA NEWTONA 1 -----------------" << endl;
        neres1 = newton_method(-10., 10., func1, dfunc1);
        cout << endl;
        cout << "----------------- METODA NEWTONA 2 -----------------" << endl;
        neres2 = newton_method(-10., 10., func2, dfunc2);
    }
    catch (const runtime_error& e) {
        cerr << "Błąd w metodzie Newtona: " << e.what() << endl;
    }

    try {
        cout << endl;
        cout << "----------------- METODA PICARDA 1 -----------------" << endl;
        pires1 = picard_method(-10., 10., func1, phi_func1, phi_func1_d);
        cout << endl;
        cout << "----------------- METODA PICARDA 2 -----------------" << endl;
        pires2 = picard_method(-10., 10., func2, phi_func2, phi_func2_d);
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
    cout << "metoda Picarda\t\t " << pires1 << "\t\t " << pires2 << endl;
    return 0;
}