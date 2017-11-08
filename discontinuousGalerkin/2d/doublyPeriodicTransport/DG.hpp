#ifndef DG_HPP
#define DG_HPP

#include <iostream>
#include <cmath>
#include "functions.hpp"

class DG {
    
    public:
        
        DG( const double&, const double&, const double&, const double&,
        const int&, const int&, const int&,
        const int&, const int&, double&, const double& );
        
        int getDFperLayer() { return n; }
        int getDF() { return N; }
        double getLayerThickness() { return dz; }
        void getXcoords( double[] );
        void getLayerMidpoints( double[] );
        
        void setGLL( double[], double[] );
        void setCardinalDerivatives( double[], double[], double[], double[] );
        
        void odeFun( double[], double[] );
        void rk( double[] );
    
    private:
        
        double a;           //left endpoint
        double b;           //right endpoint
        double c;           //bottom endpoint
        double d;           //top endpoint
        int np;             //number of polynomials per element
        int ne;             //number of elements
        int nLev;           //number of vertical levels
        int rkStages;       //number of Runge-Kutta stages (2, 3, or 4)
        int nTimesteps;     //number of time steps
        double t;           //time
        double dt;          //time increment
        
        int n;              //degrees of freedom per layer
        int N;              //total degrees of freedom
        double dz;          //layer thickness
        double* xGLL;       //GLL quadrature nodes on [-1,1]
        double* wGLL;       //GLL quadrature weights on [-1,1]
        double* xb;         //element boundaries (single layer)
        double* dx;         //element widths (single layer)
        double* xc;         //element centers (single layer)
        double* x;          //x-coordinates (single layer)
        double* z;          //array of z-coordinates (layer midpoints)
        double* weights;    //quadrature weights (single layer)
        double* dphi0dx;    //derivative of cardinal function phi0
        double* dphi1dx;    //derivative of cardinal function phi1
        double* dphi2dx;    //derivative of cardinal function phi2
        double* dphi3dx;    //derivative of cardinal function phi3
        double* u;          //horizontal velocity
        double* w;          //vertical velocity
        double* s1;         //RK stage 1
        double* s2;         //RK stage 2
        double* s3;         //RK stage 3
        double* s4;         //RK stage 4
        double* tmp;
        double alphaMax;    //Lax-Friedrichs flux parameter
        double tmpD;
        int i, j, k;
};

DG::DG( const double& A, const double& B, const double& C, const double& D,
const int& NP, const int& NE, const int& NLEV,
const int& RKSTAGES, const int& NTIMESTEPS, double& T, const double& DT ) {
    
    a = A;
    b = B;
    c = C;
    d = D;
    if( NP == 2 || NP == 3 || NP == 4 ) {
        np = NP;
    }
    else {
        std::cerr << "Error: np should be 2, 3, or 4.";
        std::exit( EXIT_FAILURE );
    }
    ne = NE;
    nLev = NLEV;
    if( RKSTAGES == 2 || RKSTAGES == 3 || RKSTAGES == 4 ) {
        rkStages = RKSTAGES;
    }
    else {
        std::cerr << "Error: rkStages should be 2, 3, or 4.";
        std::exit( EXIT_FAILURE );
    }
    nTimesteps = NTIMESTEPS;
    t = T;
    dt = DT;
    
    ///////////////////////////////////////////////////////////////////////
    
    //number of nodes (degrees of freedom) per level:
    n = np * ne;

    //total number of degrees of freedom:
    N = n * nLev;
    
    //GLL nodes and weights on [-1,1]:
    xGLL = new double[np];
    wGLL = new double[np];
    setGLL( xGLL, wGLL );
    
    //element boundaries:
    xb = new double[ne+1];
    for( i=0; i<ne+1; i++ ) {
        xb[i] = a + i*(b-a)/ne;
    }
    
    //element widths and centers:
    dx = new double[ne];
    xc = new double[ne];
    for( i=0; i<ne; i++ ) {
        dx[i] = xb[i+1] - xb[i];
        xc[i] = ( xb[i] + xb[i+1] ) / 2.;
    }
    
    //nodes and quadrature weights in a single layer:
    x = new double[n];
    weights = new double[n];
    for( i=0; i<ne; i++ ) {
        for( j=0; j<np; j++ ) {
            x[np*i+j] = xc[i] + dx[i]/2. * xGLL[j];
            weights[np*i+j] = dx[i]/2. * wGLL[j];
        }
    }
    
    //derivatives of cardinal functions:
    dphi0dx = new double[n];
    dphi1dx = new double[n];
    dphi2dx = new double[n];
    dphi3dx = new double[n];
    setCardinalDerivatives( dphi0dx, dphi1dx, dphi2dx, dphi3dx );
    
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
        for( i=0; i<nLev; i++ ) {
            z[i] = ( c + dz/2. ) + i * dz;
        }
    }
    
    //horizontal and vertical velocity:
    u = new double[N];
    w = new double[N];
    
    //arrays needed for RK time stepping:
    s1 = new double[N];
    s2 = new double[N];
    s3 = new double[N];
    s4 = new double[N];
    tmp = new double[N];
    
}

//Accessors (getters):

inline void DG::getXcoords( double X[] ) {
    for( int i=0; i<n; i++ ) {
        X[i] = x[i];
    }
}

inline void DG::getLayerMidpoints( double Z[] ) {
    for( int i=0; i<nLev; i++ ) {
        Z[i] = z[i];
    }
}

//Setters:

inline void DG::setGLL( double xGLL[], double wGLL[] ) {
    if( np == 2 ) {
        xGLL[0] = -1.;  xGLL[1] = 1.;
        wGLL[0] = 1.;  wGLL[1] = 1.;
    }
    else if( np == 3 ) {
        xGLL[0] = -1.;  xGLL[1] = 0.;  xGLL[2] = 1.;
        wGLL[0] = 1./3.;  wGLL[1] = 4./3.;  wGLL[2] = 1./3.;
    }
	else if( np == 4 ) {
        xGLL[0] = -1.;  xGLL[1] = -sqrt(5.)/5.;  xGLL[2] = sqrt(5.)/5;  xGLL[3] = 1.;
        wGLL[0] = 1./6.;  wGLL[1] = 5./6.;  wGLL[2] = 5./6.;  wGLL[3] = 1./6.;
    }
}

inline void DG::setCardinalDerivatives( double dphi0dx[], double dphi1dx[], double dphi2dx[], double dphi3dx[] ) {
    double x0, x1, x2, x3;
    for( int i=0; i<ne; i++ ) {
        if( np == 2 ) {
            x0 = x[np*i];
            x1 = x[np*i+1];
            
            dphi0dx[np*i]   = 1./(x0-x1);
            dphi0dx[np*i+1] = 1./(x0-x1);
            
            dphi1dx[np*i]   = 1./(x1-x0);
            dphi1dx[np*i+1] = 1./(x1-x0);
        }
        else if( np == 3) {
            x0 = x[np*i];
            x1 = x[np*i+1];
            x2 = x[np*i+2];
            
            dphi0dx[np*i]   = ( (x0-x1) + (x0-x2) ) / ( (x0-x1)*(x0-x2) );
            dphi0dx[np*i+1] = (x1-x2) / ( (x0-x1)*(x0-x2) );
            dphi0dx[np*i+2] = (x2-x1) / ( (x0-x1)*(x0-x2) );
            
            dphi1dx[np*i]   = (x0-x2) / ( (x1-x0)*(x1-x2) );
            dphi1dx[np*i+1] = ( (x1-x0) + (x1-x2) ) / ( (x1-x0)*(x1-x2) );
            dphi1dx[np*i+2] = (x2-x0) / ( (x1-x0)*(x1-x2) );
            
            dphi2dx[np*i]   = (x0-x1) / ( (x2-x0)*(x2-x1) );
            dphi2dx[np*i+1] = (x1-x0) / ( (x2-x0)*(x2-x1) );
            dphi2dx[np*i+2] = ( (x2-x0) + (x2-x1) ) / ( (x2-x0)*(x2-x1) );
        }
        else if( np == 4 ) {
            x0 = x[np*i];
            x1 = x[np*i+1];
            x2 = x[np*i+2];
            x3 = x[np*i+3];
            
            dphi0dx[np*i]   = ( (x0-x1)*(x0-x2) + (x0-x1)*(x0-x3) + (x0-x2)*(x0-x3) )
            / ( (x0-x1)*(x0-x2)*(x0-x3) );
            dphi0dx[np*i+1] = (x1-x2)*(x1-x3) / ( (x0-x1)*(x0-x2)*(x0-x3) );
            dphi0dx[np*i+2] = (x2-x1)*(x2-x3) / ( (x0-x1)*(x0-x2)*(x0-x3) );
            dphi0dx[np*i+3] = (x3-x1)*(x3-x2) / ( (x0-x1)*(x0-x2)*(x0-x3) );
            
            dphi1dx[np*i]   = (x0-x2)*(x0-x3) / ( (x1-x0)*(x1-x2)*(x1-x3) );
            dphi1dx[np*i+1] = ( (x1-x0)*(x1-x2) + (x1-x0)*(x1-x3) + (x1-x2)*(x1-x3) )
            / ( (x1-x0)*(x1-x2)*(x1-x3) );
            dphi1dx[np*i+2] = (x2-x0)*(x2-x3) / ( (x1-x0)*(x1-x2)*(x1-x3) );
            dphi1dx[np*i+3] = (x3-x0)*(x3-x2) / ( (x1-x0)*(x1-x2)*(x1-x3) );
            
            dphi2dx[np*i]   = (x0-x1)*(x0-x3) / ( (x2-x0)*(x2-x1)*(x2-x3) );
            dphi2dx[np*i+1] = (x1-x0)*(x1-x3) / ( (x2-x0)*(x2-x1)*(x2-x3) ); 
            dphi2dx[np*i+2] = ( (x2-x0)*(x2-x1) + (x2-x0)*(x2-x3) + (x2-x1)*(x2-x3) )
            / ( (x2-x0)*(x2-x1)*(x2-x3) );
            dphi2dx[np*i+3] = (x3-x0)*(x3-x1) / ( (x2-x0)*(x2-x1)*(x2-x3) );
            
            dphi3dx[np*i]   = (x0-x1)*(x0-x2) / ( (x3-x0)*(x3-x1)*(x3-x2) );
            dphi3dx[np*i+1] = (x1-x0)*(x1-x2) / ( (x3-x0)*(x3-x1)*(x3-x2) );
            dphi3dx[np*i+2] = (x2-x0)*(x2-x1) / ( (x3-x0)*(x3-x1)*(x3-x2) );
            dphi3dx[np*i+3] = ( (x3-x0)*(x3-x1) + (x3-x0)*(x3-x2) + (x3-x1)*(x3-x2) )
            / ( (x3-x0)*(x3-x1)*(x3-x2) );
        }
    }
}

//Functions related to time stepping:

inline void DG::odeFun( double rho[], double rhoPrime[] ) {
    //For each element, you have np ODEs:  (drho/dt)_i = RHS_i, for i = 0,...,np-1.
    //This function computes RHS.  The output is the 1D array rhoPrime.
    //Set parameter for the Lax-Friedrichs flux (LFF) for DG:
    alphaMax = 0.;
    for( i=0; i<N; i++ ) {
        tmpD = std::abs( u[i] );
        if( tmpD > alphaMax ) { alphaMax = tmpD; }
    }
    for( k=0; k<nLev; k++ ) {
        //lateral operations using DG:
        for( i=0; i<ne; i++ ) {
            if( np > 2 ) { rhoPrime[k*n+(np*i+1)] = 0; }
            if( np > 3 ) { rhoPrime[k*n+(np*i+2)] = 0; }
            if( i == 0 ) { //left-most node of left-most element (periodic BC enforcement and LFF):
                rhoPrime[k*n+(np*i)] = ( rho[k*n+(n-1)]*u[k*n+(n-1)] + rho[k*n+(np*i)]*u[k*n+(np*i)] )/2.
                - alphaMax * ( rho[k*n+(np*i)] - rho[k*n+(n-1)] );
            }
            else { //Lax-Friedrichs Flux (LFF) for the left-most node of all other elements:
                rhoPrime[k*n+(np*i)] = ( rho[k*n+(np*i-1)]*u[k*n+(np*i-1)]
                + rho[k*n+(np*i)]*u[k*n+(np*i)] )/2.
                - alphaMax * ( rho[k*n+(np*i)] - rho[k*n+(np*i-1)] );
            }
            if( i == ne-1 ) { //right-most node of right-most element (periodic BC enforcement and LFF):
                rhoPrime[k*n+(np*i+np-1)] = -( ( rho[k*n+(np*i+np-1)]*u[k*n+(np*i+np-1)]
                + rho[k*n+(0)]*u[k*n+(0)] )/2.
                - alphaMax * ( rho[k*n+(0)] - rho[k*n+(np*i+np-1)] ) );
            }
            else { //Lax-Friedrichs Flux (LFF) for the right-most node of all other elements:
                rhoPrime[k*n+(np*i+np-1)] = -( ( rho[k*n+(np*i+np-1)]*u[k*n+(np*i+np-1)]
                + rho[k*n+(np*i+np)]*u[k*n+(np*i+np)] )/2.
                - alphaMax * ( rho[k*n+(np*i+np)] - rho[k*n+(np*i+np-1)] ) );
            }
            for( j=0; j<np; j++ ) { //the non-flux part of the RHS:
                tmpD = weights[np*i+j] * ( rho[k*n+(np*i+j)] * u[k*n+(np*i+j)] ); 
                rhoPrime[k*n+(np*i)]   = rhoPrime[k*n+(np*i)]   + tmpD * dphi0dx[np*i+j];
                rhoPrime[k*n+(np*i+1)] = rhoPrime[k*n+(np*i+1)] + tmpD * dphi1dx[np*i+j];
                if( np > 2 ) { rhoPrime[k*n+(np*i+2)] = rhoPrime[k*n+(np*i+2)] + tmpD * dphi2dx[np*i+j]; }
                if( np > 3 ) { rhoPrime[k*n+(np*i+3)] = rhoPrime[k*n+(np*i+3)] + tmpD * dphi3dx[np*i+j]; }
            }
            for( j=0; j<np; j++ ) {
                rhoPrime[k*n+(np*i+j)] = rhoPrime[k*n+(np*i+j)] / weights[np*i+j];
            }
        }
        //vertical operations using FD2:
        if( nLev != 1 ) {
            if( k == 0 ) {
                for( i=0; i<n; i++ ) {
                    rhoPrime[k*n+i] = rhoPrime[k*n+i]
                    - ( rho[(k+1)*n+i]*w[(k+1)*n+i] - rho[(nLev-1)*n+i]*w[(nLev-1)*n+i] ) / (2*dz);
                }
            }
            else if( k == nLev-1 ) {
                for( i=0; i<n; i++ ) {
                    rhoPrime[k*n+i] = rhoPrime[k*n+i]
                    - ( rho[0*n+i]*w[0*n+i] - rho[(k-1)*n+i]*w[(k-1)*n+i] ) / (2*dz);
                }
            }
            else {
                for( i=0; i<n; i++ ) {
                    rhoPrime[k*n+i] = rhoPrime[k*n+i]
                    - ( rho[(k+1)*n+i]*w[(k+1)*n+i] - rho[(k-1)*n+i]*w[(k-1)*n+i] ) / (2*dz);
                }
            }
        }
    }
}

inline void DG::rk( double rho[] ) {
    if( rkStages == 2 ) {
        //2-stage, 2nd order RK.  The outputs are t and rho.  t increments and rho gets updated.
        //stage 1:
        odeFun( rho, s1 );
        //stage 2:
        t = t + dt/2.;
        for( j=0; j<nLev; j++ ) {
            for( i=0; i<n; i++ ) {
                tmp[j*n+i] = rho[j*n+i] + dt/2. * s1[j*n+i];
                u[j*n+i] = uFunc( x[i], z[j], t );
                w[j*n+i] = wFunc( x[i], z[j], t );
            }
        }
        odeFun( tmp, s1 );
        //update t and get new value of rho:
        t = t + dt/2.;
        for( i=0; i<N; i++ ) {
            rho[i] = rho[i] + dt * s1[i];
        }
    }
    else if( rkStages == 3 ) {
        //3-stage, 3rd order RK.  The outputs are t and rho.  t increments and rho gets updated.
        //stage 1:
        odeFun( rho, s1 );
        //stage 2:
        t = t + dt/3.;
        for( j=0; j<nLev; j++ ) {
            for( i=0; i<n; i++ ) {
                tmp[j*n+i] = rho[j*n+i] + dt/3. * s1[j*n+i];
                u[j*n+i] = uFunc( x[i], z[j], t );
                w[j*n+i] = wFunc( x[i], z[j], t );
            }
        }
        odeFun( tmp, s2 );
        //stage 3:
        t = t + dt/3.;
        for( j=0; j<nLev; j++ ) {
            for( i=0; i<n; i++ ) {
                tmp[j*n+i] = rho[j*n+i] + 2*dt/3. * s2[j*n+i];
                u[j*n+i] = uFunc( x[i], z[j], t );
                w[j*n+i] = wFunc( x[i], z[j], t );
            }
        }
        odeFun( tmp, s2 );
        //update t and get new value of rho:
        t = t + dt/3.;
        for( i=0; i<N; i++ ) {
            rho[i] = rho[i] + dt/4. * ( s1[i] + 3*s2[i] );
        } 
    }
    else if( rkStages == 4 ) {
        //4-stage, 4th order RK.  The outputs are t and rho.  t increments and rho gets updated.
        //stage 1:
        odeFun( rho, s1 );
        //stage 2:
        t = t + dt/2.;
        for( j=0; j<nLev; j++ ) {
            for( i=0; i<n; i++ ) {
                tmp[j*n+i] = rho[j*n+i] + dt/2. * s1[j*n+i];
                u[j*n+i] = uFunc( x[i], z[j], t );
                w[j*n+i] = wFunc( x[i], z[j], t );
            }
        }
        odeFun( tmp, s2 );
        //stage 3:
        for( i=0; i<N; i++ ) {
            tmp[i] = rho[i] + dt/2. * s2[i];
        }
        odeFun( tmp, s3 );
        //stage 4:
        t = t + dt/2;
        for( j=0; j<nLev; j++ ) {
            for( i=0; i<n; i++ ) {
                tmp[j*n+i] = rho[j*n+i] + dt * s3[j*n+i];
                u[j*n+i] = uFunc( x[i], z[j], t );
                w[j*n+i] = wFunc( x[i], z[j], t );
            }
        }
        odeFun( tmp, s4 );
        //get new value:
        for( i=0; i<N; i++ ) {
            rho[i] = rho[i] + dt/6. * ( s1[i] + 2*s2[i] + 2*s3[i] + s4[i] );
        }
    }
}

#endif
