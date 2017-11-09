#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <cmath>

//all of these should be periodic laterally.

//Function describing the initial condition for rho:
double rhoIC( const double& x, const double& z ) {
    return exp( -10 * ( pow(x,2) + pow(z,2) ) );
}

//Function describing the initial condition for rhoU:
double rhoUic( const double& x, const double& z ) {
    return 0.;
}

//Function describing the initial condition for the density rho:
double rhoWic( const double& x, const double& z ) {
    return 0.;
}

//Function describing the initial condition for the density rho:
double rhoThIC( const double& x, const double& z ) {
    return 0.;
}

#endif
