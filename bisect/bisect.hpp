#ifndef BISECT_HPP
#define BISECT_HPP

#include <iostream>
#include <iomanip>
#include "f.hpp"

double bisect( double xLeft, double xRight, const double& epsilon ) {
    double xMid = ( xLeft + xRight ) / 2.;

    if( xRight-xLeft < epsilon ) {
    }
    else if( f(xMid) == 0. ) {
    }
    else if( f(xLeft) * f(xMid) < 0. ) {
        xMid = bisect( xLeft, xMid, epsilon );
    }
    else {
        xMid = bisect( xMid, xRight, epsilon );
    }
    std::cout << "(" << std::setprecision(10) << std::setw(12) << xLeft << ", " << xRight << ")" << std::endl;
    return xMid;
}

#endif
