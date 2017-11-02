#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <iostream>
#include <cmath>

//Function describing the initial condition for rho, which should be periodic:
double rhoIC( double& x ) {
    //return exp( -10 * pow(x,2) );
    return cos( M_PI * x );
    //if( ( x > -1./2. ) && ( x < 1./2. ) ) {
    //    return 1./2.;
    //}
    //else {
    //    return 0.;
    //}
}

//Function describing the given velocity u:
double uFunc( double& x, double& t ) {
    //if( t <= 1. ) {
    //    return 1.;
    //}
    //else {
    //    return -1.;
    //}
    //return sin( M_PI * t );
    return 1.;
    //return exp( -(t-1)*(t-1) );
}

///////////////////////////////////////////////////////////////////////////

void getElementBoundaries( const int& ne, const double& a, const double& b, double xb[] ) {
    for( int i=0; i<ne+1; i++ ) {
        xb[i] = a + i*(b-a)/ne;
    }
}

void getElementWidthsAndCenters( const int& ne, double xb[], double dx[], double xc[] ) {
    for( int i=0; i<ne; i++ ) {
        dx[i] = xb[i+1] - xb[i];
        xc[i] = ( xb[i] + xb[i+1] ) / 2.;
    }
}

int getGLL( const int& np, double xGLL[], double wGLL[] ) {
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
    else {
        std::cerr << "Error:  np should be 2, 3, or 4.";
        return EXIT_FAILURE;
    }
    return 0;
}

void getCoordsAndQuadWeights( const int& ne, const int& np, double xc[], double dx[],
double xGLL[], double wGLL[], double x[], double w[] ) {
    for( int i=0; i<ne; i++ ) {
        for( int j=0; j<np; j++ ) {
            x[np*i+j] = xc[i] + dx[i]/2. * xGLL[j];
            w[np*i+j] = dx[i]/2. * wGLL[j];
        }
    }
}

int getCardinalDerivatives( const int& ne, const int& np, double x[], double dphi0dx[],
double dphi1dx[], double dphi2dx[], double dphi3dx[] ) {
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
        else {
            std::cerr << "Error:  np should be 2, 3, or 4.";
            return EXIT_FAILURE;
        }
    }
    return 0;
}

void odeFun( int& i, int& j, const int& ne, const int& np, const int& N, double u[], double x[], double& t,
double& alphaMax, double w[], double rho[], double rhoPrime[], double dphi0dx[], double dphi1dx[],
double dphi2dx[], double dphi3dx[], double& tmpD )
//For each element, you have np ODEs:  (drho/dt)_i = RHS_i, for i = 0,...,np-1.
//This function computes RHS.  The output is the 1D array rhoPrime.
{
    alphaMax = 0.;
    for( i=0; i<N; i++ ) {
        if( std::abs(u[i]) > alphaMax ) {
            alphaMax = std::abs(u[i]);
        }
    }
    for( i=0; i<ne; i++ ) {
        if( np > 2 ) { rhoPrime[np*i+1] = 0; }
        if( np > 3 ) { rhoPrime[np*i+2] = 0; }
        //Left-most element:
        if( i == 0 ) { //left-most node of left-most element (periodic BC enforcement and LFF):
            rhoPrime[np*i] = ( rho[N-1]*u[N-1] + rho[np*i]*u[np*i] )/2.
            - alphaMax * ( rho[np*i] - rho[N-1] );
        }
        else { //Lax-Friedrichs flux (LFF):
            rhoPrime[np*i] = ( rho[np*i-1]*u[np*i-1] + rho[np*i]*u[np*i] )/2.
            - alphaMax * ( rho[np*i] - rho[np*i-1] );
        }
        //Right-most element:
        if( i == ne-1 ) { //right-most node of right-most element (periodic BC enforcement and LFF):
            rhoPrime[np*i+(np-1)] = -( ( rho[np*i+(np-1)]*u[np*i+(np-1)] + rho[0]*u[0] )/2.
            - alphaMax * ( rho[0] - rho[np*i+(np-1)] ) );
        }
        else { //Lax-Friedrichs flux (LFF):
            rhoPrime[np*i+(np-1)] = -( ( rho[np*i+(np-1)]*u[np*i+(np-1)] + rho[np*i+np]*u[np*i+np] )/2.
            - alphaMax * ( rho[np*i+np] - rho[np*i+(np-1)] ) );
        }
        for( j=0; j<np; j++ ) {
            tmpD = w[np*i+j]*(rho[np*i+j]*u[np*i+j]); 
            rhoPrime[np*i]   = rhoPrime[np*i]   + tmpD * dphi0dx[np*i+j];
            rhoPrime[np*i+1] = rhoPrime[np*i+1] + tmpD * dphi1dx[np*i+j];
            if( np > 2 ) {
                rhoPrime[np*i+2] = rhoPrime[np*i+2] + tmpD * dphi2dx[np*i+j];
            }
            if( np > 3 ) {
                rhoPrime[np*i+3] = rhoPrime[np*i+3] + tmpD * dphi3dx[np*i+j];
            }
        }
        for( j=0; j<np; j++ ) {
            rhoPrime[np*i+j] = rhoPrime[np*i+j] / w[np*i+j];
        }
    }
}

void rk3( int& i, int& j, const int& ne, const int& np, const int& N, double u[], double x[], double& t,
double& alphaMax, double w[], double rho[], double dphi0dx[], double dphi1dx[], double dphi2dx[],
double dphi3dx[], const double& dt,  double s1[], double s2[], double s3[], double s4[], double tmp[],
double& tmpD ) {
    //3-stage, 3rd order RK.  The output is the 1D array rho.
    //stage 1:
    odeFun( i, j, ne, np, N, u, x, t, alphaMax, w, rho, s1, dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //stage 2:
    t = t + dt/3.;
    for( i=0; i<N; i++ ) {
        tmp[i] = rho[i] + dt/3.*s1[i];
        u[i] = uFunc( x[i], t );
    }
    odeFun( i, j, ne, np, N, u, x, t, alphaMax, w, tmp, s2, dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //stage 3:
    t = t + dt/3.;
    for( i=0; i<N; i++ ) {
        tmp[i] = rho[i] + 2*dt/3.*s2[i];
        u[i] = uFunc( x[i], t );
    }
    odeFun( i, j, ne, np, N, u, x, t, alphaMax, w, tmp, s2, dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //update t and get new value of rho:
    t = t + dt/3.;
    for( i=0; i<N; i++ ) {
        rho[i] = rho[i] + dt/4. * ( s1[i] + 3*s2[i] );
    }
}

void rk4( int& i, int& j, const int& ne, const int& np, const int& N, double u[], double x[], double& t,
double& alphaMax, double w[], double rho[], double dphi0dx[], double dphi1dx[], double dphi2dx[],
double dphi3dx[], const double& dt,  double s1[], double s2[], double s3[], double s4[], double tmp[],
double& tmpD ) {
    //4-stage, 4th order RK.  The output is the 1D array rho.
    //stage 1:
    odeFun( i, j, ne, np, N, u, x, t, alphaMax, w, rho, s1, dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //stage 2:
    t = t + dt/2.;
    for( i=0; i<N; i++ ) {
        tmp[i] = rho[i] + dt/2.*s1[i];
        u[i] = uFunc( x[i], t );
    }
    odeFun( i, j, ne, np, N, u, x, t, alphaMax, w, tmp, s2, dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //stage 3:
    for( i=0; i<N; i++ ) {
        tmp[i] = rho[i] + dt/2.*s2[i];
    }
    odeFun( i, j, ne, np, N, u, x, t, alphaMax, w, tmp, s3, dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //stage 4:
    t = t + dt/2;
    for( i=0; i<N; i++ ) {
        tmp[i] = rho[i] + dt*s3[i];
        u[i] = uFunc( x[i], t );
    }
    odeFun( i, j, ne, np, N, u, x, t, alphaMax, w, tmp, s4, dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //get new value:
    for( i=0; i<N; i++ ) {
        rho[i] = rho[i] + dt/6. * ( s1[i] + 2*s2[i] + 2*s3[i] + s4[i] );
    }
}

#endif
