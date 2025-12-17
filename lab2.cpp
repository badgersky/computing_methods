#include <iostream>
#include <fstream>
#include <cmath>
#include <filesystem>

long double func(double x) {
    double x2 = x * x;
    double x4 = x2 * x2;
    double x6 = x4 * x2;
    double x8 = x4 * x4;
    double x10 = x8 * x2;
    if (fabs(x) < 1e-1) {

        return 1. / (1. + (x2 / 20. + x4 / 840. + x6 / 60480 + x8 / 6652800 + x10 / 1037836800));
    }

    return (x * x * x) / (6. * (sinh(x) - x));
}

int main() {
    std::ifstream dataf("C:\\Users\\szymon\\workbench\\computing_methods\\dane_do_laboratorium_2.txt");
    if (!dataf.is_open()) {
        std::cerr << "lol";
        return 1;
    }
    
    std::ofstream resf("C:\\Users\\szymon\\workbench\\computing_methods\\wyniki_porownania_lab_2.txt");
    if (!resf.is_open()) {
        std::cerr << "lol";
        return 1;
    }

    std::string line;
    std::getline(dataf, line);
    std::getline(dataf, line);
    std::getline(dataf, line);

    long double logx, x_exact, err, log_err;
    double x;
    while(dataf >> logx >> x >> x_exact) {
        err = (func(x) - x_exact) / x_exact;
        log_err = log10l(fabsl(err));
        
        resf << logx << " " << log_err << std::endl;
    }

    dataf.close();
    return 0;
}