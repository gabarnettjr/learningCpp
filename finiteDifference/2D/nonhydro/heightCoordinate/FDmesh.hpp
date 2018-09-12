#ifndef FDMESH_HPP
#define FDMESH_HPP

#include <iostream>
#include <cmath>

#include "Parameters.hpp"

class FDmesh {
    
    public:
        
        FDmesh( const std::string&, const Parameters& );
        
        double* x;
        double* s;
        
    private:
        
        double* zSurfB;
        double* zSurfPrimeB;
        double* dzB;
        double* zSurf;       //surface
        double* zSurfPrime;  //derivative of surface
        double* dz;          //layer thicknesses
        double* xb;          //x coordinates of cell boundaries
        double* zb;          //z coordinates of cell boundaries
        double* z;           //all z coordinates
};

FDmesh::FDmesh( const std::string& st, const Parameters& P ) {
    
    //array of x coordinates of cell boundaries on a single level:
    xb = new double[P.n+1];
    for( int i = 0; i < P.n+1; i++ ) {
        xb[i] = P.a + i * P.dx;
    }
    
    // averaging to get array of x-coordinates in one layer:
    x = new double[P.n];
    for( int i = 0; i < P.n; i++ ) {
        x[i] = ( xb[i] + xb[i+1] ) / 2.;
    }
    
    //stuff related to topography:
    zSurfB = new double[P.n+1];
    zSurfPrimeB = new double[P.n+1];
    dzB = new double[P.n+1];
    for( int i = 0; i < P.n+1; i++ ) {
        if( st == "risingBubble" ) {
            zSurfB[i] = 500 * ( 1. + sin(2*M_PI*xb[i]/5000.) )
            zSurfPrimeB[i] = M_PI/5. * cos(2*M_PI*xb[i]/5000.);
            dzB[i] = ( P.d - zSurfB[i] ) / P.nLev;
        }
    }
    
    //averaging:
    zSurf = new double[P.n];
    zSurfPrime = new double[P.n];
    dz = new double[P.n];
    for( int i = 0; i < P.n; i++ ) {
        zSurf[i] = ( zSurfB[i] + zSurfB[i+1] ) / 2.;
        zSurfPrime[i] = ( zSurfPrimeB[i] + zSurfPrimeB[i+1] ) / 2.;
        dz[i] = ( dzB[i] + dzB[i+1] ) / 2.;
    }
    
    //array of all layer midpoints (z coordinates):
    z = new double[P.N];
    for( int j = 0; j < P.nLev; j++ ) {
        for( int i = 0; i < P.n; i++ ) {
            z[j*P.n+i] = ( zSurf[i] + P.dz[i]/2. ) + j * P.dz[i];
        }
    }
}

#endif
