#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <cmath>

//Function describing the initial condition for rho, which should be periodic:
double rhoIC( const double& x, const double& z ) {
    return exp( -10 * ( pow(x,2) + pow(z,2) ) );
    // return cos( M_PI * x ) + sin( M_PI * z );
}

//Function describing the prescribed velocity u:
double uFunc( const double& x, const double& z, const double& t ) {
    // return 1.;
    return sin( M_PI * t );
}

double wFunc( const double& x, const double& z, const double& t ) {
    return 1.;
    // return cos( M_PI * t );
}

#endif
