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
        void getXcoords( double[] );
        void getLayerMidpoints( double[] );
        
        void setGLL( double[], double[] );
        void setCardinalDerivatives( double[], double[], double[], double[] );
        
        void singleVariableRHS( const double[], const double[], const double[], double[] );
        void odeFun( const double[], const double[], const double[], const double[],
        double[], double[], double[], double[] );
        void rk( double[], double[], double[], double[] );
    
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
        double* P;          //Pressure (diagnostic)
        double* s1_rho;     //RK stage 1
        double* s1_rhoU;
        double* s1_rhoW;
        double* s1_rhoTh;
        double* s2_rho;     //RK stage 2
        double* s2_rhoU;
        double* s2_rhoW;
        double* s2_rhoTh;
        double* s3_rho;     //RK stage 3
        double* s3_rhoU;
        double* s3_rhoW;
        double* s3_rhoTh;
        double* s4_rho;     //RK stage 4
        double* s4_rhoU;
        double* s4_rhoW;
        double* s4_rhoTh;
        double* tmp_rho;    //used inside rk
        double* tmp_rhoU;
        double* tmp_rhoW;
        double* tmp_rhoTh;
        double alphaMax;    //Lax-Friedrichs flux parameter
        double tmpD;        //used inside odeFun
        double* F;
        double* G;
        int i, j, k;        //used in for-loops
        double Po;             //physical constants
        double Cp;
        double Cv;
        double Rd;
        double g;
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
    
    //arrays needed for RK time stepping:
    s1_rho    = new double[N];
    s1_rhoU   = new double[N];
    s1_rhoW   = new double[N];
    s1_rhoTh  = new double[N];
    s2_rho    = new double[N];
    s2_rhoU   = new double[N];
    s2_rhoW   = new double[N];
    s2_rhoTh  = new double[N];
    s3_rho    = new double[N];
    s3_rhoU   = new double[N];
    s3_rhoW   = new double[N];
    s3_rhoTh  = new double[N];
    s4_rho    = new double[N];
    s4_rhoU   = new double[N];
    s4_rhoW   = new double[N];
    s4_rhoTh  = new double[N];
    tmp_rho   = new double[N];
    tmp_rhoU  = new double[N];
    tmp_rhoW  = new double[N];
    tmp_rhoTh = new double[N];
    
    //physical constants:
    Po = 100000.;
    Cp = 1004.;
    Cv = 717.;
    Rd = Cp - Cv;
    g = 9.81;
    
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

inline void DG::singleVariableRHS( const double rho[], const double rhoU[], const double rhoW[],
double rhoPrime[] ) {
    for( k=0; k<nLev; k++ ) {
        //lateral operations using DG:
        for( i=0; i<ne; i++ ) {
            if( np > 2 ) { rhoPrime[k*n+(np*i+1)] = 0; }
            if( np > 3 ) { rhoPrime[k*n+(np*i+2)] = 0; }
            if( i == 0 ) { //left-most node of left-most element (periodic BC enforcement and LFF):
                rhoPrime[k*n+(np*i)] = ( rhoU[k*n+(n-1)] + rhoU[k*n+(np*i)] )/2.
                - alphaMax * ( rho[k*n+(np*i)] - rho[k*n+(n-1)] );
            }
            else { //Lax-Friedrichs Flux (LFF) for the left-most node of all other elements:
                rhoPrime[k*n+(np*i)] = ( rhoU[k*n+(np*i-1)] + rhoU[k*n+(np*i)] )/2.
                - alphaMax * ( rho[k*n+(np*i)] - rho[k*n+(np*i-1)] );
            }
            if( i == ne-1 ) { //right-most node of right-most element (periodic BC enforcement and LFF):
                rhoPrime[k*n+(np*i+np-1)] = -( ( rhoU[k*n+(np*i+np-1)] + rhoU[k*n+(0)] )/2.
                - alphaMax * ( rho[k*n+(0)] - rho[k*n+(np*i+np-1)] ) );
            }
            else { //Lax-Friedrichs Flux (LFF) for the right-most node of all other elements:
                rhoPrime[k*n+(np*i+np-1)] = -( ( rhoU[k*n+(np*i+np-1)] + rhoU[k*n+(np*i+np)] )/2.
                - alphaMax * ( rho[k*n+(np*i+np)] - rho[k*n+(np*i+np-1)] ) );
            }
            for( j=0; j<np; j++ ) { //the non-flux part of the RHS:
                tmpD = weights[np*i+j] * rhoU[k*n+(np*i+j)]; 
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
        if( k == 0 ) {
            for( i=0; i<n; i++ ) { //reflect over bottom boundary since w=0 there:
                rhoPrime[k*n+i] = rhoPrime[k*n+i] - ( 2 * rhoW[(k+1)*n+i] ) / (2*dz);
            }
        }
        else if( k == nLev-1 ) { //reflect over top boundary since w=0 there:
            for( i=0; i<n; i++ ) {
                rhoPrime[k*n+i] = rhoPrime[k*n+i] - ( -2 * rhoW[(k-1)*n+i] ) / (2*dz);
            }
        }
        else {
            for( i=0; i<n; i++ ) {
                rhoPrime[k*n+i] = rhoPrime[k*n+i] - ( rhoW[(k+1)*n+i] - rhoW[(k-1)*n+i] ) / (2*dz);
            }
        }
    }
}

inline void DG::odeFun( const double rho[], const double rhoU[], const double rhoW[], const double rhoTh[],
double rhoPrime[], double rhoUprime[], double rhoWprime[], double rhoThPrime[] ) {
    alphaMax = 0.;
    for( i=0; i<N; i++ ) {
        tmpD = std::abs( rhoU[i] / rho[i] );
        if( tmpD > alphaMax ) { alphaMax = tmpD; }
        P[i] = Po * pow( Rd*rhoTh[i]/Po, Cp/Cv );
    }
    //rho:
    singleVariableRHS( rho, rhoU, rhoW, rhoPrime );
    //rhoU:
    for( i=0; i<N; i++ ) {
        F[i] = rhoU[i] * rhoU[i] / rho[i] + P[i];
        G[i] = rhoU[i] * rhoW[i] / rho[i];
    }
    singleVariableRHS( rhoU, F, G, rhoUprime );
    //rhoW:
    for( i=0; i<N; i++ ) {
        F[i] = rhoW[i] * rhoW[i] / rho[i];
    }
    singleVariableRHS( rhoW, G, F, rhoWprime );
    for( k=0; k<nLev; k++ ) { //adding missing terms ( - dP/dz - rho*g ):
        if( k == 0 ) {
            for( i=0; i<n; i++ ) { //lowest layer:
                tmpD = P[k*n+i] + g*dz/2 * ( 3*rho[k*n+i] - rho[(k+1)*n+i] );
                rhoWprime[k*n+i] = rhoWprime[k*n+i] - ( P[(k+1)*n+i] - tmpD ) / (2*dz)
                - rho[k*n+i] * g;
            }
        }
        else if( k == nLev-1 ) { //highest layer:
            for( i=0; i<n; i++ ) {
                tmpD = P[k*n+i] - g*dz/2 * ( 3*rho[k*n+i] - rho[(k-1)*n+i] );
                rhoWprime[k*n+i] = rhoWprime[k*n+i] - ( tmpD - P[(k-1)*n+i] ) / (2*dz)
                - rho[k*n+i] * g;
            }
        }
        else {
            for( i=0; i<n; i++ ) {
                rhoWprime[k*n+i] = rhoWprime[k*n+i] - ( P[(k+1)*n+i] - P[(k-1)*n+i] ) / (2*dz)
                - rho[k*n+i] * g;
            }
        }
    }
    //rhoTh:
    for( i=0; i<N; i++ ) {
        F[i] = rhoU[i] * rhoTh[i] / rho[i];
        G[i] = rhoW[i] * rhoTh[i] / rho[i];
    }
    singleVariableRHS( rhoTh, F, G, rhoThPrime );
}

inline void DG::rk( double rho[], double rhoU[], double rhoW[], double rhoTh[] ) {
    if( rkStages == 2 ) {
        //2-stage, 2nd order RK.  The outputs are t and rho.  t increments and rho gets updated.
        //stage 1:
        odeFun( rho, rhoU, rhoW, rhoTh, s1_rho, s1_rhoU, s1_rhoW, s1_rhoTh );
        //stage 2:
        t = t + dt/2.;
        for( j=0; j<nLev; j++ ) {
            for( i=0; i<n; i++ ) {
                tmp_rho[j*n+i]   = rho[j*n+i]   + dt/2. * s1_rho[j*n+i];
                tmp_rhoU[j*n+i]  = rhoU[j*n+i]  + dt/2. * s1_rhoU[j*n+i];
                tmp_rhoW[j*n+i]  = rhoW[j*n+i]  + dt/2. * s1_rhoW[j*n+i];
                tmp_rhoTh[j*n+i] = rhoTh[j*n+i] + dt/2. * s1_rhoTh[j*n+i];
            }
        }
        odeFun( tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s1_rho, s1_rhoU, s1_rhoW, s1_rhoTh );
        //update t and get new value of rho:
        t = t + dt/2.;
        for( i=0; i<N; i++ ) {
            rho[i]   = rho[i]   + dt * s1_rho[i];
            rhoU[i]  = rhoU[i]  + dt * s1_rhoU[i];
            rhoW[i]  = rhoW[i]  + dt * s1_rhoW[i];
            rhoTh[i] = rhoTh[i] + dt * s1_rhoTh[i];
        }
    }
    else if( rkStages == 3 ) {
        //3-stage, 3rd order RK.  The outputs are t and rho.  t increments and rho gets updated.
        //stage 1:
        odeFun( rho, rhoU, rhoW, rhoTh, s1_rho, s1_rhoU, s1_rhoW, s1_rhoTh );
        //stage 2:
        t = t + dt/3.;
        for( j=0; j<nLev; j++ ) {
            for( i=0; i<n; i++ ) {
                tmp_rho[j*n+i]   = rho[j*n+i]   + dt/3. * s1_rho[j*n+i];
                tmp_rhoU[j*n+i]  = rhoU[j*n+i]  + dt/3. * s1_rhoU[j*n+i];
                tmp_rhoW[j*n+i]  = rhoW[j*n+i]  + dt/3. * s1_rhoW[j*n+i];
                tmp_rhoTh[j*n+i] = rhoTh[j*n+i] + dt/3. * s1_rhoTh[j*n+i];
            }
        }
        odeFun( tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s2_rho, s2_rhoU, s2_rhoW, s2_rhoTh );
        //stage 3:
        t = t + dt/3.;
        for( j=0; j<nLev; j++ ) {
            for( i=0; i<n; i++ ) {
                tmp_rho[j*n+i]   = rho[j*n+i]   + 2*dt/3. * s2_rho[j*n+i];
                tmp_rhoU[j*n+i]  = rhoU[j*n+i]  + 2*dt/3. * s2_rhoU[j*n+i];
                tmp_rhoW[j*n+i]  = rhoW[j*n+i]  + 2*dt/3. * s2_rhoW[j*n+i];
                tmp_rhoTh[j*n+i] = rhoTh[j*n+i] + 2*dt/3. * s2_rhoTh[j*n+i];
            }
        }
        odeFun( tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s2_rho, s2_rhoU, s2_rhoW, s2_rhoTh );
        //update t and get new value of rho:
        t = t + dt/3.;
        for( i=0; i<N; i++ ) {
            rho[i]   = rho[i]   + dt/4. * ( s1_rho[i]   + 3*s2_rho[i] );
            rhoU[i]  = rhoU[i]  + dt/4. * ( s1_rhoU[i]  + 3*s2_rhoU[i] );
            rhoW[i]  = rhoW[i]  + dt/4. * ( s1_rhoW[i]  + 3*s2_rhoW[i] );
            rhoTh[i] = rhoTh[i] + dt/4. * ( s1_rhoTh[i] + 3*s2_rhoTh[i] );
        } 
    }
    else if( rkStages == 4 ) {
        //4-stage, 4th order RK.  The outputs are t and rho.  t increments and rho gets updated.
        //stage 1:
        odeFun( rho, rhoU, rhoW, rhoTh, s1_rho, s1_rhoU, s1_rhoW, s1_rhoTh );
        //stage 2:
        t = t + dt/2.;
        for( j=0; j<nLev; j++ ) {
            for( i=0; i<n; i++ ) {
                tmp_rho[j*n+i]   = rho[j*n+i]   + dt/2. * s1_rho[j*n+i];
                tmp_rhoU[j*n+i]  = rhoU[j*n+i]  + dt/2. * s1_rhoU[j*n+i];
                tmp_rhoW[j*n+i]  = rhoW[j*n+i]  + dt/2. * s1_rhoW[j*n+i];
                tmp_rhoTh[j*n+i] = rhoTh[j*n+i] + dt/2. * s1_rhoTh[j*n+i];
            }
        }
        odeFun( tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s2_rho, s2_rhoU, s2_rhoW, s2_rhoTh );
        //stage 3:
        for( i=0; i<N; i++ ) {
            tmp_rho[i]   = rho[i]   + dt/2. * s2_rho[i];
            tmp_rhoU[i]  = rhoU[i]  + dt/2. * s2_rhoU[i];
            tmp_rhoW[i]  = rhoW[i]  + dt/2. * s2_rhoW[i];
            tmp_rhoTh[i] = rhoTh[i] + dt/2. * s2_rhoTh[i];
        }
        odeFun( tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s3_rho, s3_rhoU, s3_rhoW, s3_rhoTh );
        //stage 4:
        t = t + dt/2;
        for( j=0; j<nLev; j++ ) {
            for( i=0; i<n; i++ ) {
                tmp_rho[j*n+i]   = rho[j*n+i]   + dt * s3_rho[j*n+i];
                tmp_rhoU[j*n+i]  = rhoU[j*n+i]  + dt * s3_rhoU[j*n+i];
                tmp_rhoW[j*n+i]  = rhoW[j*n+i]  + dt * s3_rhoW[j*n+i];
                tmp_rhoTh[j*n+i] = rhoTh[j*n+i] + dt * s3_rhoTh[j*n+i];
            }
        }
        odeFun( tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s4_rho, s4_rhoU, s4_rhoW, s4_rhoTh );
        //get new value:
        for( i=0; i<N; i++ ) {
            rho[i]   = rho[i]   + dt/6. * ( s1_rho[i]   + 2*s2_rho[i]   + 2*s3_rho[i]   + s4_rho[i] );
            rhoU[i]  = rhoU[i]  + dt/6. * ( s1_rhoU[i]  + 2*s2_rhoU[i]  + 2*s3_rhoU[i]  + s4_rhoU[i] );
            rhoW[i]  = rhoW[i]  + dt/6. * ( s1_rhoW[i]  + 2*s2_rhoW[i]  + 2*s3_rhoW[i]  + s4_rhoW[i] );
            rhoTh[i] = rhoTh[i] + dt/6. * ( s1_rhoTh[i] + 2*s2_rhoTh[i] + 2*s3_rhoTh[i] + s4_rhoTh[i] );
        }
    }
}

#endif
