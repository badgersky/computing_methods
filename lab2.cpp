#include <iostream>
#include <fstream>
#include <cmath>
#include <filesystem>

long double func(double x) {
    if (fabsl(x) < 1e-2) {
        return 1. / (1. + (x * x / 20. + x * x * x * x / 840. + x * x * x * x * x * x / 60480));
    } 

    return (x * x * x) / (6. * (sinhl(x) - x));
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

    long double logx, x, x_exact, err, log_err;
    while(dataf >> logx >> x >> x_exact) {
        err = (func(x) - x_exact) / x_exact;
        log_err = log10l(fabsl(err));
        
        resf << logx << " " << log_err << std::endl;
    }

    dataf.close();
    return 0;
}