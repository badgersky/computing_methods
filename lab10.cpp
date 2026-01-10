#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

#define N 10000

using namespace std;

double f(double t, double y) {
    return  -((100 * t + 10) / (t + 1)) * (y - 1);
}

double BME(double t, double y, double h) {
    return (y + h * f(t, y));
}

double PME(double t, double y, double h) {
    double t2 = t + h;
    return ((y + h * ((100 * t2 + 10) / (t2 + 1))) / (1 + h * ((100 * t2 + 10) / (t2 + 1))));
}

double PMT(double t, double y, double h) {
    double t2 = t + h;
    return (((y + h * f(t, y) / 2) + (h * (100 * t2 + 10) / (2 * (t2 + 1)))) / (1 + h * (100 * t2 + 10) / (2 * (t2 + 1))));
}

double exact(double t) {
    return 1 + pow(1 + t, 90) * exp(-100. * t);
}

int main() {
    ofstream f1("C:\\Users\\szymon\\workbench\\computing_methods\\results_lab10.txt");
    ofstream f2("C:\\Users\\szymon\\workbench\\computing_methods\\errors_lab10.txt");
    ofstream f3("C:\\Users\\szymon\\workbench\\computing_methods\\bme_unstable_lab10.txt");

    for (int n = 10; n <= N; n++) {
        double t_s = 0.;
        double t_e = 1.;
        double y = 2.;
        double y_bme, y_pme, y_pmt, y_e;
        double h = fabs(t_s - t_e) / (n - 1);
        double t = t_s;
        double y1 = y, y2 = y, y3 = y;
        double err_bme = 0.;
        double err_pme = 0.;
        double err_pmt = 0.;

        if (!f1.is_open() || !f2.is_open()) {
            return 1;
        }

        for (int i = 0; i < n; i++) {
            t = t_s + i * h;
            y_e = exact(t);
            
            if (n == 1000) {
                f1 << setw(18) << setprecision(15) << fixed << y1 << "\t" << y2 << "\t" << y3 << "\t" << y_e << endl;
            }

            if (n == 10) {
                f3 << setw(18) << setprecision(15) << fixed << y1 << "\t" << y_e << endl;
            }

            y_bme = BME(t, y1, h);
            y_pme = PME(t, y2, h);
            y_pmt = PMT(t, y3, h);
            
            if (abs(y1 - y_e) > err_bme) err_bme = abs(y1 - y_e);
            if (abs(y2 - y_e) > err_pme) err_pme = abs(y2 - y_e);
            if (abs(y3 - y_e) > err_pmt) err_pmt = abs(y3 - y_e); 

            y1 = y_bme; y2 = y_pme, y3 = y_pmt;
        }

        f2 << setw(18) << setprecision(15) << fixed << h << "\t" << err_bme << "\t" << err_pme << "\t" << err_pmt << "\t" << endl;
    }

    f1.close();
    f2.close();
    f3.close();
}