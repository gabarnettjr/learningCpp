#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <iostream>
#include <cmath>

int getCardinalDerivatives( const int& ne, const int& np, double x[], double dphi0dx[], double dphi1dx[], double dphi2dx[], double dphi3dx[] ) {
    double x0;
    double x1;
    double x2;
    double x3;
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
            
            dphi0dx[np*i]   = ( (x0-x1)*(x0-x2) + (x0-x1)*(x0-x3) + (x0-x2)*(x0-x3) ) / ( (x0-x1)*(x0-x2)*(x0-x3) );
            dphi0dx[np*i+1] = (x1-x2)*(x1-x3) / ( (x0-x1)*(x0-x2)*(x0-x3) );
            dphi0dx[np*i+2] = (x2-x1)*(x2-x3) / ( (x0-x1)*(x0-x2)*(x0-x3) );
            dphi0dx[np*i+3] = (x3-x1)*(x3-x2) / ( (x0-x1)*(x0-x2)*(x0-x3) );
            
            dphi1dx[np*i]   = (x0-x2)*(x0-x3) / ( (x1-x0)*(x1-x2)*(x1-x3) );
            dphi1dx[np*i+1] = ( (x1-x0)*(x1-x2) + (x1-x0)*(x1-x3) + (x1-x2)*(x1-x3) ) / ( (x1-x0)*(x1-x2)*(x1-x3) );
            dphi1dx[np*i+2] = (x2-x0)*(x2-x3) / ( (x1-x0)*(x1-x2)*(x1-x3) );
            dphi1dx[np*i+3] = (x3-x0)*(x3-x2) / ( (x1-x0)*(x1-x2)*(x1-x3) );
            
            dphi2dx[np*i]   = (x0-x1)*(x0-x3) / ( (x2-x0)*(x2-x1)*(x2-x3) );
            dphi2dx[np*i+1] = (x1-x0)*(x1-x3) / ( (x2-x0)*(x2-x1)*(x2-x3) ); 
            dphi2dx[np*i+2] = ( (x2-x0)*(x2-x1) + (x2-x0)*(x2-x3) + (x2-x1)*(x2-x3) ) / ( (x2-x0)*(x2-x1)*(x2-x3) );
            dphi2dx[np*i+3] = (x3-x0)*(x3-x1) / ( (x2-x0)*(x2-x1)*(x2-x3) );
            
            dphi3dx[np*i]   = (x0-x1)*(x0-x2) / ( (x3-x0)*(x3-x1)*(x3-x2) );
            dphi3dx[np*i+1] = (x1-x0)*(x1-x2) / ( (x3-x0)*(x3-x1)*(x3-x2) );
            dphi3dx[np*i+2] = (x2-x0)*(x2-x1) / ( (x3-x0)*(x3-x1)*(x3-x2) );
            dphi3dx[np*i+3] = ( (x3-x0)*(x3-x1) + (x3-x0)*(x3-x2) + (x3-x1)*(x3-x2) ) / ( (x3-x0)*(x3-x1)*(x3-x2) );
        }
        else {
            std::cerr << "Error:  np should be 2, 3, or 4.";
            return EXIT_FAILURE;
            // throw std::logic_error( "np should be 2, 3, or 4." );
        }
    }
    return 0;
}

void odeFun( int i, int j, const int& ne, const int& np, const int& N, const double& u, double t, double w[], double rho[], double rhoPrime[], double dphi0dx[], double dphi1dx[], double dphi2dx[], double dphi3dx[] )
//For each element, you have np ODEs:  (drho/dt)_i = RHS_i, for i = 0,...,np-1.
//This function computes RHS.  The output is the 1D array rhoPrime.
{
    for( i=0; i<ne; i++ ) {
        rhoPrime[np*i] = 0;
        rhoPrime[np*i+1] = 0;
        if( np > 2 ) { rhoPrime[np*i+2] = 0; }
        if( np > 3 ) { rhoPrime[np*i+3] = 0; }
        //Left-most element:
        if( i == 0 ) { //left-most node of left-most element (periodic BC enforcement and LFF):
            rhoPrime[np*i] = u*( rho[N-1] + rho[np*i] )/2. - std::abs(u) * ( rho[np*i] - rho[N-1] );
        }
        else { //Lax-Friedrichs flux (LFF):
            rhoPrime[np*i] = u*( rho[np*i-1] + rho[np*i] )/2. - std::abs(u) * ( rho[np*i] - rho[np*i-1] );
        }
        //Right-most element:
        if( i == ne-1 ) { //right-most node of right-most element (periodic BC enforcement and LFF):
            rhoPrime[np*i+(np-1)] = -( u*( rho[np*i+(np-1)] + rho[0] )/2. - std::abs(u) * ( rho[0] - rho[np*i+(np-1)] ) );
        }
        else { //Lax-Friedrichs flux (LFF):
            rhoPrime[np*i+(np-1)] = -( u*( rho[np*i+(np-1)] + rho[np*i+np] )/2. - std::abs(u) * ( rho[np*i+np] - rho[np*i+(np-1)] ) );
        }
        for( j=0; j<np; j++ ) {
            t = w[np*i+j] * (rho[np*i+j]*u);            //temporary value (time t is not explicitly needed)
            rhoPrime[np*i]   = rhoPrime[np*i]   + t * dphi0dx[np*i+j];
            rhoPrime[np*i+1] = rhoPrime[np*i+1] + t * dphi1dx[np*i+j];
            if( np > 2 ) { rhoPrime[np*i+2] = rhoPrime[np*i+2] + t * dphi2dx[np*i+j]; }
            if( np > 3 ) { rhoPrime[np*i+3] = rhoPrime[np*i+3] + t * dphi3dx[np*i+j]; }
        }
        for( j=0; j<np; j++ ) {
            rhoPrime[np*i+j] = rhoPrime[np*i+j] / w[np*i+j];
        }
    }
}

void rk( int i, int j, const int& ne, const int& np, const int& N, const double& u, double t, double w[], double rho[], double dphi0dx[], double dphi1dx[], double dphi2dx[], double dphi3dx[], const double& dt, double s1[], double s2[], double s3[], double s4[], double tmp[] ) {
    //The output is the 1D array rho.  Currently only does RK4.
    odeFun( i, j, ne, np, N, u, t,       w, rho, s1, dphi0dx, dphi1dx, dphi2dx, dphi3dx );
    for( i=0; i<N; i++ ) {
        tmp[i] = rho[i] + dt/2.*s1[i];
    }
    odeFun( i, j, ne, np, N, u, t+dt/2., w, tmp, s2, dphi0dx, dphi1dx, dphi2dx, dphi3dx );
    for( i=0; i<N; i++ ) {
        tmp[i] = rho[i] + dt/2.*s2[i];
    }
    odeFun( i, j, ne, np, N, u, t+dt/2., w, tmp, s3, dphi0dx, dphi1dx, dphi2dx, dphi3dx );
    for( i=0; i<N; i++ ) {
        tmp[i] = rho[i] + dt*s3[i];
    }
    odeFun( i, j, ne, np, N, u, t+dt,    w, tmp, s4, dphi0dx, dphi1dx, dphi2dx, dphi3dx );
    for( i=0; i<N; i++ ) {
        rho[i] = rho[i] + dt/6. * ( s1[i] + 2*s2[i] + 2*s3[i] + s4[i] );
    }
}

#endif