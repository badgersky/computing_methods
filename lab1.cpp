#include <iostream>
#include <cmath>

int main() {
    float x1 = 1.0f;
    int ft1 = 0;
    while (1.0f + x1 > 1.0f) {
        x1 /= 2.0f;
        ft1++;
    }
    ft1--;
  
    std::cout << "liczba bitów mantysy dla float:                              " << ft1 << std::endl;
    std::cout << "liczba cyfr znaczących dla float:                            " << std::ceil((ft1 * std::log10(2))) << std::endl;
    std::cout << "znaleziony epsilon maszynowy dla float:                      " << (x1*2.0f) << std::endl;
    std::cout << "epsilon maszynowy dla float z biblioteki standardowej:       " << std::numeric_limits<float>::epsilon() << std::endl;

    double x2 = 1.0;
    int ft2 = 0;
    while (1.0 + x2 > 1.0) {
        x2 /= 2.0;
        ft2++;
    }
    ft2--;

    std::cout << "------------------------------------------------------------------------" << std::endl;
    std::cout << "liczba bitów mantysy dla double:                             " << ft2 << std::endl;
    std::cout << "liczba cyfr znaczących dla double:                           " << std::ceil((ft2 * std::log10(2))) << std::endl;
    std::cout << "znaleziony epsilon maszynowy dla double:                     " << (x2*2.0) << std::endl;
    std::cout << "epsilon maszynowy dla double z biblioteki standardowej:      " << std::numeric_limits<double>::epsilon() << std::endl;

    long double x3 = 1.0l;
    int ft3 = 0;
    while (1.0l + x3 > 1.0l) {
        x3 /= 2.0l;
        ft3++;
    }
    ft3--;

    std::cout << "------------------------------------------------------------------------" << std::endl;
    std::cout << "liczba bitów mantysy dla long double:                        " << ft3 << std::endl;
    std::cout << "liczba cyfr znaczących dla long double:                      " << std::ceil((ft3 * std::log10(2))) << std::endl;
    std::cout << "znaleziony epsilon maszynowy dla long double:                " << (x3*2.0l) << std::endl;
    std::cout << "epsilon maszynowy dla long double z biblioteki standardowej: " << std::numeric_limits<long double>::epsilon() << std::endl;
}