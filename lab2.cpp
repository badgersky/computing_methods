#include <iostream>
#include <fstream>
#include <cmath>
#include <filesystem>

double func(double x) {
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

    double logx, x, x_exact, err, log_err;
    while(dataf >> logx >> x >> x_exact) {
        err = fabs((func(x) - x_exact) / x_exact);
        log_err = log10(fabs(err));
        
        resf << logx << " " << log_err << std::endl;
    }

    dataf.close();
    return 0;
}