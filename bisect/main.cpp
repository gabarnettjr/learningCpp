#include <iostream>
#include <iomanip>
#include <cmath>
#include "bisect.hpp"

int main()
{
    const double xLeft = 0.;
    const double xRight = 10.;
    const double epsilon = .0000001;
    
    double root = bisect( xLeft, xRight, epsilon );
    double exact = sqrt(2.);

    int pr = 10;

    std::cout << "The approximation is:" << std::endl << std::setprecision(pr) << root << std::endl;
    std::cout << std::endl;
    std::cout << "The actual root is:" << std::endl << std::setprecision(pr) << exact << std::endl;
    std::cout << std::endl;
    std::cout << "The absolute error is:" << std::endl << std::setprecision(pr) << std::abs(root-exact) << std::endl;

    return 0;
}
