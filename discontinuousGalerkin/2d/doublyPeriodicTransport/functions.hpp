#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <iostream>
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

///////////////////////////////////////////////////////////////////////////

void odeFun( int& i, int& j, int& k, const int& ne, const int& np, const int& nLev, const int& n,
const int& N, const double u[], const double w[], const double x[], const double z[], const double& t,
const double& dz, double& alphaMax, const double weights[], const double rho[], double rhoPrime[],
const double dphi0dx[], const double dphi1dx[], const double dphi2dx[], const double dphi3dx[],
double& tmpD ) {
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

void rk2( int& i, int& j, int& k, const int& ne, const int& np, const int& nLev, const int& n, const int& N,
double u[], double w[], const double x[], const double z[], double& t, const double& dz,
double& alphaMax, const double weights[], double rho[],
const double dphi0dx[], const double dphi1dx[], const double dphi2dx[], const double dphi3dx[],
const double& dt,  double s1[], double s2[], double s3[], double s4[], double tmp[], double& tmpD ) {
    //2-stage, 2nd order RK.  The outputs are t and rho.  t increments and rho gets updated.
    //stage 1:
    odeFun( i, j, k, ne, np, nLev, n, N, u, w, x, z, t, dz, alphaMax, weights, rho, s1,
    dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //stage 2:
    t = t + dt/2.;
    for( j=0; j<nLev; j++ ) {
        for( i=0; i<n; i++ ) {
            tmp[j*n+i] = rho[j*n+i] + dt/2. * s1[j*n+i];
            u[j*n+i] = uFunc( x[i], z[j], t );
            w[j*n+i] = wFunc( x[i], z[j], t );
        }
    }
    odeFun( i, j, k, ne, np, nLev, n, N, u, w, x, z, t, dz, alphaMax, weights, tmp, s1,
    dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //update t and get new value of rho:
    t = t + dt/2.;
    for( i=0; i<N; i++ ) {
        rho[i] = rho[i] + dt * s1[i];
    }
}

void rk3( int& i, int& j, int& k, const int& ne, const int& np, const int& nLev, const int& n, const int& N,
double u[], double w[], const double x[], const double z[], double& t, const double& dz,
double& alphaMax, const double weights[], double rho[],
const double dphi0dx[], const double dphi1dx[], const double dphi2dx[], const double dphi3dx[],
const double& dt,  double s1[], double s2[], double s3[], double s4[], double tmp[], double& tmpD ) {
    //3-stage, 3rd order RK.  The outputs are t and rho.  t increments and rho gets updated.
    //stage 1:
    odeFun( i, j, k, ne, np, nLev, n, N, u, w, x, z, t, dz, alphaMax, weights, rho, s1,
    dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //stage 2:
    t = t + dt/3.;
    for( j=0; j<nLev; j++ ) {
        for( i=0; i<n; i++ ) {
            tmp[j*n+i] = rho[j*n+i] + dt/3. * s1[j*n+i];
            u[j*n+i] = uFunc( x[i], z[j], t );
            w[j*n+i] = wFunc( x[i], z[j], t );
        }
    }
    odeFun( i, j, k, ne, np, nLev, n, N, u, w, x, z, t, dz, alphaMax, weights, tmp, s2,
    dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //stage 3:
    t = t + dt/3.;
    for( j=0; j<nLev; j++ ) {
        for( i=0; i<n; i++ ) {
            tmp[j*n+i] = rho[j*n+i] + 2*dt/3. * s2[j*n+i];
            u[j*n+i] = uFunc( x[i], z[j], t );
            w[j*n+i] = wFunc( x[i], z[j], t );
        }
    }
    odeFun( i, j, k, ne, np, nLev, n, N, u, w, x, z, t, dz, alphaMax, weights, tmp, s2,
    dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //update t and get new value of rho:
    t = t + dt/3.;
    for( i=0; i<N; i++ ) {
        rho[i] = rho[i] + dt/4. * ( s1[i] + 3*s2[i] );
    }
}

void rk4( int& i, int& j, int& k, const int& ne, const int& np, const int& nLev, const int& n, const int& N,
double u[], double w[], const double x[], const double z[], double& t, const double& dz,
double& alphaMax, const double weights[], double rho[],
const double dphi0dx[], const double dphi1dx[], const double dphi2dx[], const double dphi3dx[],
const double& dt,  double s1[], double s2[], double s3[], double s4[], double tmp[], double& tmpD ) {
    //4-stage, 4th order RK.  The outputs are t and rho.  t increments and rho gets updated.
    //stage 1:
    odeFun( i, j, k, ne, np, nLev, n, N, u, w, x, z, t, dz, alphaMax, weights, rho, s1,
    dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //stage 2:
    t = t + dt/2.;
    for( j=0; j<nLev; j++ ) {
        for( i=0; i<n; i++ ) {
            tmp[j*n+i] = rho[j*n+i] + dt/2. * s1[j*n+i];
            u[j*n+i] = uFunc( x[i], z[j], t );
            w[j*n+i] = wFunc( x[i], z[j], t );
        }
    }
    odeFun( i, j, k, ne, np, nLev, n, N, u, w, x, z, t, dz, alphaMax, weights, tmp, s2,
    dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //stage 3:
    for( i=0; i<N; i++ ) {
        tmp[i] = rho[i] + dt/2. * s2[i];
    }
    odeFun( i, j, k, ne, np, nLev, n, N, u, w, x, z, t, dz, alphaMax, weights, tmp, s3,
    dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //stage 4:
    t = t + dt/2;
    for( j=0; j<nLev; j++ ) {
        for( i=0; i<n; i++ ) {
            tmp[j*n+i] = rho[j*n+i] + dt * s3[j*n+i];
            u[j*n+i] = uFunc( x[i], z[j], t );
            w[j*n+i] = wFunc( x[i], z[j], t );
        }
    }
    odeFun( i, j, k, ne, np, nLev, n, N, u, w, x, z, t, dz, alphaMax, weights, tmp, s4,
    dphi0dx, dphi1dx, dphi2dx, dphi3dx, tmpD );
    //get new value:
    for( i=0; i<N; i++ ) {
        rho[i] = rho[i] + dt/6. * ( s1[i] + 2*s2[i] + 2*s3[i] + s4[i] );
    }
}

#endif
