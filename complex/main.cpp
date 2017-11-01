#include <iostream>
#include "complex.hpp"

int main()
{
    double a = 1.;
    double b = -3.;
    Complex c( a, b );

    double alpha = 9.;
    double beta = 8.;
    Complex d( alpha, beta );
    
    Complex e = c * d;
    
    std::cout << "re: " << e.getRe() << std::endl;
    std::cout << "im: " << e.getIm() << std::endl;

    e -= c;

    std::cout << "re: " << e.getRe() << std::endl;
    std::cout << "im: " << e.getIm() << std::endl;
    
    return 0;
}
