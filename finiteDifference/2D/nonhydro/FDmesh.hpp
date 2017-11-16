#ifndef FDMESH_HPP
#define FDMESH_HPP

#include <iostream>
#include <cmath>

#include "Parameters.hpp"

class FDmesh {
    
    public:
        
        FDmesh( const Parameters& );
        
        double* x;          //x coordinates in a single layer
        double* z;          //z coordinates in a single column 
};

FDmesh::FDmesh( const Parameters& P ) {
    
    //array of cell centers (x coordinates):
    x  = new double[P.n];
    for( int i = 0; i < P.n; i++ ) {
        x[i]  = ( P.a + P.dx/2 ) + i * P.dx;
    }
    
    //array of layer midpoints (z coordinates):
    z = new double[P.nLev];
    for( int i = 0; i < P.nLev; i++ ) {
        z[i] = ( P.c + P.dz/2. ) + i * P.dz;
    }
}

#endif
