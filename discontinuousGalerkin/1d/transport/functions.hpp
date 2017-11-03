#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <iostream>
#include <cmath>

//Function describing the initial condition for rho, which should be periodic:
double rhoIC( const double& x ) {
    return exp( -10 * pow(x,2) );
    //return cos( M_PI * x );
    //if( ( x > -1./2. ) && ( x < 1./2. ) ) {
    //    return 1./2.;
    //}
    //else {
    //    return 0.;
    //}
}

//Function describing the given velocity u:
double uFunc( const double& x, const double& t ) {
    // if( t < 1. ) {
        // return 1.;
    // }
    // else {
        // return -1.;
    // }
    return sin( M_PI * t );
    //return 1.;
    //return exp( -(t-1)*(t-1) );
}

///////////////////////////////////////////////////////////////////////////

void odeFun( int& i, int& j, const int& ne, const int& np, const int& N, const double u[], const double x[],
const double& t, double& alphaMax, const double w[], const double rho[], double rhoPrime[], const double dphi0dx[],
const double dphi1dx[], const double dphi2dx[], const double dphi3dx[], double& tmpD ) {
    //For each element, you have np ODEs:  (drho/dt)_i = RHS_i, for i = 0,...,np-1.
    //This function computes RHS.  The output is the 1D array rhoPrime.
    alphaMax = 0.;
    for( i=0; i<N; i++ ) {
        tmpD = std::abs( u[i] );
        if( tmpD > alphaMax ) { alphaMax = tmpD; }
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
            if( np > 2 ) { rhoPrime[np*i+2] = rhoPrime[np*i+2] + tmpD * dphi2dx[np*i+j]; }
            if( np > 3 ) { rhoPrime[np*i+3] = rhoPrime[np*i+3] + tmpD * dphi3dx[np*i+j]; }
        }
        for( j=0; j<np; j++ ) {
            rhoPrime[np*i+j] = rhoPrime[np*i+j] / w[np*i+j];
        }
    }
}

void rk2( int& i, int& j, const int& ne, const int& np, const int& N, double u[], const double x[], double& t,
double& alphaMax, const double w[], double rho[], const double dphi0dx[], const double dphi1dx[], const double dphi2dx[],
const double dphi3dx[], const double& dt,  double s1[], double s2[], double s3[], double s4[], double tmp[],
double& tmpD ) {
    //2-stage, 2nd order RK.  The outputs are t and rho.  t increments and rho gets updated.
    //stage 1:
    odeFun( i, j, ne, np, N, u, x, t, alphaMax, w, rho, s1, dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //stage 2:
    t = t + dt/2.;
    for( i=0; i<N; i++ ) {
        tmp[i] = rho[i] + dt/2.*s1[i];
        u[i] = uFunc( x[i], t );
    }
    odeFun( i, j, ne, np, N, u, x, t, alphaMax, w, tmp, s1, dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //update t and get new value of rho:
    t = t + dt/2.;
    for( i=0; i<N; i++ ) {
        rho[i] = rho[i] + dt * s1[i];
    }
}

void rk3( int& i, int& j, const int& ne, const int& np, const int& N, double u[], const double x[], double& t,
double& alphaMax, const double w[], double rho[], const double dphi0dx[], const double dphi1dx[], const double dphi2dx[],
const double dphi3dx[], const double& dt,  double s1[], double s2[], double s3[], double s4[], double tmp[],
double& tmpD ) {
    //3-stage, 3rd order RK.  The outputs are t and rho.  t increments and rho gets updated.
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

void rk4( int& i, int& j, const int& ne, const int& np, const int& N, double u[], const double x[], double& t,
double& alphaMax, const double w[], double rho[], const double dphi0dx[], const double dphi1dx[], const double dphi2dx[],
const double dphi3dx[], const double& dt,  double s1[], double s2[], double s3[], double s4[], double tmp[],
double& tmpD ) {
    //4-stage, 4th order RK.  The outputs are t and rho.  t increments and rho gets updated.
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
