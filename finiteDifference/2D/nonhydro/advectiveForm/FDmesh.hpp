#ifndef FDMESH_HPP
#define FDMESH_HPP

#include <iostream>
#include <cmath>

class FDmesh {
    
    public:
        
        FDmesh( const double&, const double&, const double&, const double&,
        const int&, const int& );
        
        double a;           //left endpoint
        double b;           //right endpoint
        double c;           //bottom endpoint
        double d;           //top endpoint
        
        int nLev;           //number of vertical levels
        
        int n;              //degrees of freedom (nodes) per layer
        int N;              //total degrees of freedom (nodes)
        
        double* x;          //x coordinates in a single layer
        double* z;          //z coordinates in a single column
        
        double dx;          //cell width
        double dz;          //layer thickness
    
    private:
        
        double* xb;         //element boundaries (single layer)
        
};

FDmesh::FDmesh( const double& aa, const double& bb, const double& cc, const double& dd,
const int& nn, const int& NLEV ) {
    
    a = aa;
    b = bb;
    c = cc;
    d = dd;
    n = nn;
    nLev = NLEV;
    
    ///////////////////////////////////////////////////////////////////////

    //total number of degrees of freedom:
    N = n * nLev;
    
    //cell boundaries:
    xb = new double[n+1];
    for( int i = 0; i < n+1; i++ ) {
        xb[i] = a + i*(b-a)/n;
    }
    
    //cell width:
    dx = xb[1] - xb[0];
    
    //cell centers (x coordinates):
    x  = new double[n];
    for( int i = 0; i < n; i++ ) {
        x[i]  = ( xb[i] + xb[i+1] ) / 2.;
    }
    
    //space between layer midpoints:
    if( nLev == 1 ) {
        dz = 0.;
    }
    else {
        dz = ( d - c ) / ( nLev );
    }
    
    //array of layer midpoints:
    z = new double[nLev];
    if( nLev == 1 ) {
        z[0] = 0.;
    }
    else {
        for( int i = 0; i < nLev; i++ ) {
            z[i] = ( c + dz/2. ) + i * dz;
        }
    }
    
}

#endif
