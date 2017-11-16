#ifndef TIMESTEPPER_HPP
#define TIMESTEPPER_HPP

#include <iostream>
#include <cmath>

#include "Parameters.hpp"
#include "Variables.hpp"

class TimeStepper {
    
    public:
        
        TimeStepper( const Parameters&, Variables& );
        
        double t;
        
        void rk( const Parameters&, Variables& );
        
    private:
        
        double* p;          //temporary pressure for inside odefun
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
        double tmpD;        //used inside odeFun
        double* F;
        double* G;
        int i, j, k;        //used in for-loops
        
        void singleVariableRHS( const Parameters&, const Variables&,
        const double[], const double[], const double[],
        double[] );
        
        void odeFun( const Parameters&, const Variables&,
        const double[], const double[], const double[], const double[],
        double[], double[], double[], double[] );
};

TimeStepper::TimeStepper( const Parameters& P, Variables& V ) {
    t = P.t;
    
    //temporary pressure variable:
    p = new double[P.N];
    
    //arrays needed for RK time stepping:
    s1_rho    = new double[P.N];
    s1_rhoU   = new double[P.N];
    s1_rhoW   = new double[P.N];
    s1_rhoTh  = new double[P.N];
    s2_rho    = new double[P.N];
    s2_rhoU   = new double[P.N];
    s2_rhoW   = new double[P.N];
    s2_rhoTh  = new double[P.N];
    s3_rho    = new double[P.N];
    s3_rhoU   = new double[P.N];
    s3_rhoW   = new double[P.N];
    s3_rhoTh  = new double[P.N];
    s4_rho    = new double[P.N];
    s4_rhoU   = new double[P.N];
    s4_rhoW   = new double[P.N];
    s4_rhoTh  = new double[P.N];
    tmp_rho   = new double[P.N];
    tmp_rhoU  = new double[P.N];
    tmp_rhoW  = new double[P.N];
    tmp_rhoTh = new double[P.N];
    
    //used inside singleVariableRHS:
    F = new double[P.N];
    G = new double[P.N];
}

inline void TimeStepper::singleVariableRHS( const Parameters& P, const Variables& V,
const double rho[], const double rhoU[], const double rhoW[],
double rhoPrime[] ) {
    //Sets rhoPrime equal to -d(rhoU)/dx - d(rhoW)/dz, using FD for d/dx and d/dz.
    //This is used for each of the four prognostic variables: rho, rhoU, rhoW, rhoTh.
    for( k = 0; k < P.nLev; k++ ) {
        //lateral operations using FD:
        if( P.stenX == 3 ) {
            for( i=0; i<P.n; i++ ) {
                if( i == 0 ) {
                    rhoPrime[k*P.n+i] = - ( rhoU[k*P.n+(i+1)] - rhoU[k*P.n+(P.n-1)] ) / (2*P.dx)
                    + std::abs( V.rhoU[k*P.n+i] / V.rho[k*P.n+i] ) * P.dx/2
                    * ( rho[k*P.n+(P.n-1)] - 2*rho[k*P.n+i] + rho[k*P.n+(i+1)] ) / pow(P.dx,2.);
                }
                else if( i == P.n-1 ) {
                    rhoPrime[k*P.n+i] = - ( rhoU[k*P.n+0] - rhoU[k*P.n+(i-1)] ) / (2*P.dx)
                    + std::abs( V.rhoU[k*P.n+i] / V.rho[k*P.n+i] ) * P.dx/2
                    * ( rho[k*P.n+(i-1)] - 2*rho[k*P.n+i] + rho[k*P.n+0] ) / pow(P.dx,2.);
                }
                else {
                    rhoPrime[k*P.n+i] = - ( rhoU[k*P.n+(i+1)] - rhoU[k*P.n+(i-1)] ) / (2*P.dx)
                    + std::abs( V.rhoU[k*P.n+i] / V.rho[k*P.n+i] ) * P.dx/2
                    * ( rho[k*P.n+(i-1)] - 2*rho[k*P.n+i] + rho[k*P.n+(i+1)] ) / pow(P.dx,2.);
                }
            }
        }
        else if( P.stenX == 5 ) {
            for( i=0; i<P.n; i++ ) {
                if( i == 0 ) {
                    rhoPrime[k*P.n+i] = - ( 1./12.*rhoU[k*P.n+(P.n-2)] - 2./3.*rhoU[k*P.n+(P.n-1)]
                    + 2./3.*rhoU[k*P.n+(i+1)] - 1./12.*rhoU[k*P.n+(i+2)] ) / P.dx
                    - std::abs( V.rhoU[k*P.n+i] / V.rho[k*P.n+i] ) * pow(P.dx,3.)/12.
                    * ( rho[k*P.n+(P.n-2)] - 4.*rho[k*P.n+(P.n-1)] + 6.*rho[k*P.n+i]
                    - 4.*rho[k*P.n+(i+1)] + rho[k*P.n+(i+2)] ) / pow(P.dx,4.);
                }
                else if( i == 1 ) {
                    rhoPrime[k*P.n+i] = - ( 1./12.*rhoU[k*P.n+(P.n-1)] - 2./3.*rhoU[k*P.n+(i-1)]
                    + 2./3.*rhoU[k*P.n+(i+1)] - 1./12.*rhoU[k*P.n+(i+2)] ) / P.dx
                    - std::abs( V.rhoU[k*P.n+i] / V.rho[k*P.n+i] ) * pow(P.dx,3.)/12.
                    * ( rho[k*P.n+(P.n-1)] - 4.*rho[k*P.n+(i-1)] + 6.*rho[k*P.n+i]
                    - 4.*rho[k*P.n+(i+1)] + rho[k*P.n+(i+2)] ) / pow(P.dx,4.);
                }
                else if( i == P.n-2 ) {
                    rhoPrime[k*P.n+i] = - ( 1./12.*rhoU[k*P.n+(i-2)] - 2./3.*rhoU[k*P.n+(i-1)]
                    + 2./3.*rhoU[k*P.n+(i+1)] - 1./12.*rhoU[k*P.n+(0)] ) / P.dx
                    - std::abs( V.rhoU[k*P.n+i] / V.rho[k*P.n+i] ) * pow(P.dx,3.)/12.
                    * ( rho[k*P.n+(i-2)] - 4.*rho[k*P.n+(i-1)] + 6.*rho[k*P.n+i]
                    - 4.*rho[k*P.n+(i+1)] + rho[k*P.n+(0)] ) / pow(P.dx,4.);
                }
                else if( i == P.n-1 ) {
                    rhoPrime[k*P.n+i] = - ( 1./12.*rhoU[k*P.n+(i-2)] - 2./3.*rhoU[k*P.n+(i-1)]
                    + 2./3.*rhoU[k*P.n+(0)] - 1./12.*rhoU[k*P.n+(1)] ) / P.dx
                    - std::abs( V.rhoU[k*P.n+i] / V.rho[k*P.n+i] ) * pow(P.dx,3.)/12.
                    * ( rho[k*P.n+(i-2)] - 4.*rho[k*P.n+(i-1)] + 6.*rho[k*P.n+i]
                    - 4.*rho[k*P.n+(0)] + rho[k*P.n+(1)] ) / pow(P.dx,4.);
                }
                else {
                    rhoPrime[k*P.n+i] = - ( 1./12.*rhoU[k*P.n+(i-2)] - 2./3.*rhoU[k*P.n+(i-1)]
                    + 2./3.*rhoU[k*P.n+(i+1)] - 1./12.*rhoU[k*P.n+(i+2)] ) / P.dx
                    - std::abs( V.rhoU[k*P.n+i] / V.rho[k*P.n+i] ) * pow(P.dx,3.)/12.
                    * ( rho[k*P.n+(i-2)] - 4.*rho[k*P.n+(i-1)] + 6.*rho[k*P.n+i]
                    - 4.*rho[k*P.n+(i+1)] + rho[k*P.n+(i+2)] ) / pow(P.dx,4.);
                }
            }
        }
        else {
            std::cerr << "Error: stenX should be 3 or 5.";
            std::exit( EXIT_FAILURE );
        }
        //vertical operations using FD2:
        if( k == 0 ) {
            for( i = 0; i < P.n; i++ ) { //reflect over bottom boundary since w=0 there:
                rhoPrime[k*P.n+i] = rhoPrime[k*P.n+i] - ( rhoW[(k+1)*P.n+i] - -rhoW[k*P.n+i] ) / (2*P.dz)
                + std::abs( V.rhoW[k*P.n+i] / V.rho[k*P.n+i] ) * P.dz/2.
                * ( rho[k*P.n+i] - 2.*rho[(k+1)*P.n+i] + rho[(k+2)*P.n+i] ) / pow(P.dz,2.);
            }
        }
        else if( k == P.nLev-1 ) { //reflect over top boundary since w=0 there:
            for( i = 0; i < P.n; i++ ) {
                rhoPrime[k*P.n+i] = rhoPrime[k*P.n+i] - ( -rhoW[k*P.n+i] - rhoW[(k-1)*P.n+i] ) / (2*P.dz)
                + std::abs( V.rhoW[k*P.n+i] / V.rho[k*P.n+i] ) * P.dz/2.
                * ( rho[(k-2)*P.n+i] - 2.*rho[(k-1)*P.n+i] + rho[k*P.n+i] ) / pow(P.dz,2.);
            }
        }
        else {
            for( i = 0; i < P.n; i++ ) {
                rhoPrime[k*P.n+i] = rhoPrime[k*P.n+i] - ( rhoW[(k+1)*P.n+i] - rhoW[(k-1)*P.n+i] ) / (2*P.dz)
                + std::abs( V.rhoW[k*P.n+i] / V.rho[k*P.n+i] ) * P.dz/2.
                * ( rho[(k-1)*P.n+i] - 2.*rho[k*P.n+i] + rho[(k+1)*P.n+i] ) / pow(P.dz,2.);
            }
        }
    }
}

inline void TimeStepper::odeFun( const Parameters& P, const Variables& V,
const double rho[], const double rhoU[], const double rhoW[], const double rhoTh[],
double rhoPrime[], double rhoUprime[], double rhoWprime[], double rhoThPrime[] ) {
    for( i = 0; i < P.N; i++ ) {
        p[i] = P.Po * pow( P.Rd * rhoTh[i] / P.Po, P.Cp / P.Cv );
    }
    //rho:
    singleVariableRHS( P, V, rho, rhoU, rhoW, rhoPrime );
    //rhoU:
    for( i = 0; i < P.N; i++ ) {
        F[i] = rhoU[i] * rhoU[i] / rho[i] + p[i];
        G[i] = rhoU[i] * rhoW[i] / rho[i];
    }
    singleVariableRHS( P, V, rhoU, F, G, rhoUprime );
    //rhoW:
    for( i = 0; i < P.N; i++ ) {
        F[i] = rhoW[i] * rhoW[i] / rho[i];
    }
    singleVariableRHS( P, V, rhoW, G, F, rhoWprime );
    for( k = 0; k < P.nLev; k++ ) { //adding missing terms ( - dP/dz - rho*g ):
        if( k == 0 ) { //lowest layer (extrapolate rho and enforce dP/dz = -rho*g):
            for( i = 0; i < P.n; i++ ) {
                tmpD = p[k*P.n+i] + P.g*P.dz * ( 3.*rho[k*P.n+i] - rho[(k+1)*P.n+i] ) / 2.;
                rhoWprime[k*P.n+i] = rhoWprime[k*P.n+i] - ( p[(k+1)*P.n+i] - tmpD ) / (2.*P.dz)
                - rho[k*P.n+i] * P.g;
            }
        }
        else if( k == P.nLev-1 ) { //highest layer (extrapolate rho and enforce dP/dz = -rho*g):
            for( i = 0; i < P.n; i++ ) {
                tmpD = p[k*P.n+i] - P.g*P.dz * ( 3.*rho[k*P.n+i] - rho[(k-1)*P.n+i] ) / 2.;
                rhoWprime[k*P.n+i] = rhoWprime[k*P.n+i] - ( tmpD - p[(k-1)*P.n+i] ) / (2.*P.dz)
                - rho[k*P.n+i] * P.g;
            }
        }
        else { //all other layers:
            for( i = 0; i < P.n; i++ ) {
                rhoWprime[k*P.n+i] = rhoWprime[k*P.n+i] - ( p[(k+1)*P.n+i] - p[(k-1)*P.n+i] ) / (2.*P.dz)
                - rho[k*P.n+i] * P.g;
            }
        }
    }
    //rhoTh:
    for( i = 0; i < P.N; i++ ) {
        F[i] = rhoU[i] * rhoTh[i] / rho[i];
        G[i] = rhoW[i] * rhoTh[i] / rho[i];
    }
    singleVariableRHS( P, V, rhoTh, F, G, rhoThPrime );
}

inline void TimeStepper::rk( const Parameters& P, Variables& V ) {
    if( P.rkStages == 2 ) {
        //2-stage, 2nd order RK.  The outputs are t and rho.  t increments and rho gets updated.
        //stage 1:
        odeFun( P, V, V.rho, V.rhoU, V.rhoW, V.rhoTh, s1_rho, s1_rhoU, s1_rhoW, s1_rhoTh );
        //stage 2:
        t = t + P.dt/2.;
        for( i = 0; i < P.N; i++ ) {
            tmp_rho[i]   = V.rho[i]   + P.dt/2. * s1_rho[i];
            tmp_rhoU[i]  = V.rhoU[i]  + P.dt/2. * s1_rhoU[i];
            tmp_rhoW[i]  = V.rhoW[i]  + P.dt/2. * s1_rhoW[i];
            tmp_rhoTh[i] = V.rhoTh[i] + P.dt/2. * s1_rhoTh[i];
        }
        odeFun( P, V, tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s1_rho, s1_rhoU, s1_rhoW, s1_rhoTh );
        //update t and get new value of rho:
        t = t + P.dt/2.;
        for( i = 0; i < P.N; i++ ) {
            V.rho[i]   = V.rho[i]   + P.dt * s1_rho[i];
            V.rhoU[i]  = V.rhoU[i]  + P.dt * s1_rhoU[i];
            V.rhoW[i]  = V.rhoW[i]  + P.dt * s1_rhoW[i];
            V.rhoTh[i] = V.rhoTh[i] + P.dt * s1_rhoTh[i];
            V.p[i] = P.Po * pow( P.Rd*V.rhoTh[i]/P.Po, P.Cp/P.Cv );
        }
    }
    else if( P.rkStages == 3 ) {
        //3-stage, 3rd order RK.  The outputs are t and rho.  t increments and rho gets updated.
        //stage 1:
        odeFun( P, V, V.rho, V.rhoU, V.rhoW, V.rhoTh, s1_rho, s1_rhoU, s1_rhoW, s1_rhoTh );
        //stage 2:
        t = t + P.dt/3.;
        for( i = 0; i < P.N; i++ ) {
            tmp_rho[i]   = V.rho[i]   + P.dt/3. * s1_rho[i];
            tmp_rhoU[i]  = V.rhoU[i]  + P.dt/3. * s1_rhoU[i];
            tmp_rhoW[i]  = V.rhoW[i]  + P.dt/3. * s1_rhoW[i];
            tmp_rhoTh[i] = V.rhoTh[i] + P.dt/3. * s1_rhoTh[i];
        }
        odeFun( P, V, tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s2_rho, s2_rhoU, s2_rhoW, s2_rhoTh );
        //stage 3:
        t = t + P.dt/3.;
        for( i = 0; i < P.N; i++ ) {
            tmp_rho[i]   = V.rho[i]   + 2*P.dt/3. * s2_rho[i];
            tmp_rhoU[i]  = V.rhoU[i]  + 2*P.dt/3. * s2_rhoU[i];
            tmp_rhoW[i]  = V.rhoW[i]  + 2*P.dt/3. * s2_rhoW[i];
            tmp_rhoTh[i] = V.rhoTh[i] + 2*P.dt/3. * s2_rhoTh[i];
        }
        odeFun( P, V, tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s2_rho, s2_rhoU, s2_rhoW, s2_rhoTh );
        //update t and get new value of rho:
        t = t + P.dt/3.;
        for( i = 0; i < P.N; i++ ) {
            V.rho[i]   = V.rho[i]   + P.dt/4. * ( s1_rho[i]   + 3*s2_rho[i] );
            V.rhoU[i]  = V.rhoU[i]  + P.dt/4. * ( s1_rhoU[i]  + 3*s2_rhoU[i] );
            V.rhoW[i]  = V.rhoW[i]  + P.dt/4. * ( s1_rhoW[i]  + 3*s2_rhoW[i] );
            V.rhoTh[i] = V.rhoTh[i] + P.dt/4. * ( s1_rhoTh[i] + 3*s2_rhoTh[i] );
            V.p[i] = P.Po * pow( P.Rd*V.rhoTh[i]/P.Po, P.Cp/P.Cv );
        } 
    }
    else if( P.rkStages == 4 ) {
        //4-stage, 4th order RK.  The outputs are t and rho.  t increments and rho gets updated.
        //stage 1:
        odeFun( P, V, V.rho, V.rhoU, V.rhoW, V.rhoTh, s1_rho, s1_rhoU, s1_rhoW, s1_rhoTh );
        //stage 2:
        t = t + P.dt/2.;
        for( i = 0; i < P.N; i++ ) {
            tmp_rho[i]   = V.rho[i]   + P.dt/2. * s1_rho[i];
            tmp_rhoU[i]  = V.rhoU[i]  + P.dt/2. * s1_rhoU[i];
            tmp_rhoW[i]  = V.rhoW[i]  + P.dt/2. * s1_rhoW[i];
            tmp_rhoTh[i] = V.rhoTh[i] + P.dt/2. * s1_rhoTh[i];
        }
        odeFun( P, V, tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s2_rho, s2_rhoU, s2_rhoW, s2_rhoTh );
        //stage 3:
        for( i = 0; i < P.N; i++ ) {
            tmp_rho[i]   = V.rho[i]   + P.dt/2. * s2_rho[i];
            tmp_rhoU[i]  = V.rhoU[i]  + P.dt/2. * s2_rhoU[i];
            tmp_rhoW[i]  = V.rhoW[i]  + P.dt/2. * s2_rhoW[i];
            tmp_rhoTh[i] = V.rhoTh[i] + P.dt/2. * s2_rhoTh[i];
        }
        odeFun( P, V, tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s3_rho, s3_rhoU, s3_rhoW, s3_rhoTh );
        //stage 4:
        t = t + P.dt/2;
        for( i = 0; i < P.N; i++ ) {
            tmp_rho[i]   = V.rho[i]   + P.dt * s3_rho[i];
            tmp_rhoU[i]  = V.rhoU[i]  + P.dt * s3_rhoU[i];
            tmp_rhoW[i]  = V.rhoW[i]  + P.dt * s3_rhoW[i];
            tmp_rhoTh[i] = V.rhoTh[i] + P.dt * s3_rhoTh[i];
        }
        odeFun( P, V, tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s4_rho, s4_rhoU, s4_rhoW, s4_rhoTh );
        //get new value:
        for( i = 0; i < P.N; i++ ) {
            V.rho[i]   = V.rho[i]   + P.dt/6. * ( s1_rho[i]   + 2*s2_rho[i]   + 2*s3_rho[i]   + s4_rho[i] );
            V.rhoU[i]  = V.rhoU[i]  + P.dt/6. * ( s1_rhoU[i]  + 2*s2_rhoU[i]  + 2*s3_rhoU[i]  + s4_rhoU[i] );
            V.rhoW[i]  = V.rhoW[i]  + P.dt/6. * ( s1_rhoW[i]  + 2*s2_rhoW[i]  + 2*s3_rhoW[i]  + s4_rhoW[i] );
            V.rhoTh[i] = V.rhoTh[i] + P.dt/6. * ( s1_rhoTh[i] + 2*s2_rhoTh[i] + 2*s3_rhoTh[i] + s4_rhoTh[i] );
            V.p[i] = P.Po * pow( P.Rd*V.rhoTh[i]/P.Po, P.Cp/P.Cv );
        }
    }
    else {
        std::cerr << "Error: rkStages should be 2, 3, or 4.";
        std::exit( EXIT_FAILURE );
    }
}

#endif