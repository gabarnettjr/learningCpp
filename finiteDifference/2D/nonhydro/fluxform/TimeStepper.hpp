#ifndef TIMESTEPPER_HPP
#define TIMESTEPPER_HPP

#include <iostream>
#include <cmath>

#include "Constants.hpp"
#include "FDmesh.hpp"
#include "Variables.hpp"

class TimeStepper {
    
    public:
        
        TimeStepper( const Constants&, const FDmesh&, Variables&,
        const int&, const int&, double&, const double& );
        
        int rkStages;
        int nTimesteps;
        double t;
        double dt;
        
        void rk( const Constants&, const FDmesh&, Variables& );
        
    private:
        
        double* P;          //temporary pressure for inside odefun
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
        
        void singleVariableRHS( const Constants&, const FDmesh&, const Variables&,
        const double[], const double[], const double[],
        double[] );
        
        void odeFun( const Constants&, const FDmesh&, const Variables&,
        const double[], const double[], const double[], const double[],
        double[], double[], double[], double[] );
};

TimeStepper::TimeStepper( const Constants& C, const FDmesh& M, Variables& V,
const int& RKSTAGES, const int& NTIMESTEPS, double& T, const double& DT ) {
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
    
    //temporary pressure variable:
    P = new double[M.N];
    
    //arrays needed for RK time stepping:
    s1_rho    = new double[M.N];
    s1_rhoU   = new double[M.N];
    s1_rhoW   = new double[M.N];
    s1_rhoTh  = new double[M.N];
    s2_rho    = new double[M.N];
    s2_rhoU   = new double[M.N];
    s2_rhoW   = new double[M.N];
    s2_rhoTh  = new double[M.N];
    s3_rho    = new double[M.N];
    s3_rhoU   = new double[M.N];
    s3_rhoW   = new double[M.N];
    s3_rhoTh  = new double[M.N];
    s4_rho    = new double[M.N];
    s4_rhoU   = new double[M.N];
    s4_rhoW   = new double[M.N];
    s4_rhoTh  = new double[M.N];
    tmp_rho   = new double[M.N];
    tmp_rhoU  = new double[M.N];
    tmp_rhoW  = new double[M.N];
    tmp_rhoTh = new double[M.N];
    
    //used inside singleVariableRHS:
    F = new double[M.N];
    G = new double[M.N];
}

inline void TimeStepper::singleVariableRHS( const Constants& C, const FDmesh& M, const Variables& V,
const double rho[], const double rhoU[], const double rhoW[],
double rhoPrime[] ) {
    //Sets rhoPrime equal to -d(rhoU)/dx - d(rhoW)/dz, using FD for d/dx and d/dz.
    //This is used for each of the four prognostic variables: rho, rhoU, rhoW, rhoTh.
    for( k = 0; k < M.nLev; k++ ) {
        //lateral operations using FD:
        for( i=0; i<M.n; i++ ) {
            if( i == 0 ) {
                rhoPrime[k*M.n+i] = - ( rhoU[k*M.n+(i+1)] - rhoU[k*M.n+(M.n-1)] ) / (2*M.dx)
                + std::abs( V.rhoU[k*M.n+i] / V.rho[k*M.n+i] ) * M.dx/2
                * ( rho[k*M.n+(M.n-1)] - 2*rho[k*M.n+i] + rho[k*M.n+(i+1)] ) / pow(M.dx,2.);
            }
            else if( i == M.n-1 ) {
                rhoPrime[k*M.n+i] = - ( rhoU[k*M.n+0] - rhoU[k*M.n+(i-1)] ) / (2*M.dx)
                + std::abs( V.rhoU[k*M.n+i] / V.rho[k*M.n+i] ) * M.dx/2
                * ( rho[k*M.n+(i-1)] - 2*rho[k*M.n+i] + rho[k*M.n+0] ) / pow(M.dx,2.);
            }
            else {
                rhoPrime[k*M.n+i] = - ( rhoU[k*M.n+(i+1)] - rhoU[k*M.n+(i-1)] ) / (2*M.dx)
                + std::abs( V.rhoU[k*M.n+i] / V.rho[k*M.n+i] ) * M.dx/2
                * ( rho[k*M.n+(i-1)] - 2*rho[k*M.n+i] + rho[k*M.n+(i+1)] ) / pow(M.dx,2.);
            }
        }
        //vertical operations using FD:
        if( k == 0 ) {
            for( i = 0; i < M.n; i++ ) { //reflect over bottom boundary since w=0 there:
                rhoPrime[k*M.n+i] = rhoPrime[k*M.n+i] - ( rhoW[(k+1)*M.n+i] - -rhoW[k*M.n+i] ) / (2*M.dz);
                // + std::abs( V.rhoW[k*M.n+i] / V.rho[k*M.n+i] ) * M.dz/2.
                // * ( rho[k*M.n+i] - 2.*rho[(k+1)*M.n+i] + rho[(k+2)*M.n+i] ) / pow(M.dz,2.);
            }
        }
        else if( k == M.nLev-1 ) { //reflect over top boundary since w=0 there:
            for( i = 0; i < M.n; i++ ) {
                rhoPrime[k*M.n+i] = rhoPrime[k*M.n+i] - ( -rhoW[k*M.n+i] - rhoW[(k-1)*M.n+i] ) / (2*M.dz);
                // + std::abs( V.rhoW[k*M.n+i] / V.rho[k*M.n+i] ) * M.dz/2.
                // * ( rho[(k-2)*M.n+i] - 2.*rho[(k-1)*M.n+i] + rho[k*M.n+i] ) / pow(M.dz,2.);
            }
        }
        else {
            for( i = 0; i < M.n; i++ ) {
                rhoPrime[k*M.n+i] = rhoPrime[k*M.n+i] - ( rhoW[(k+1)*M.n+i] - rhoW[(k-1)*M.n+i] ) / (2*M.dz)
                + std::abs( V.rhoW[k*M.n+i] / V.rho[k*M.n+i] ) * M.dz/2.
                * ( rho[(k-1)*M.n+i] - 2.*rho[k*M.n+i] + rho[(k+1)*M.n+i] ) / pow(M.dz,2.);
            }
        }
    }
}

inline void TimeStepper::odeFun( const Constants& C, const FDmesh& M, const Variables& V,
const double rho[], const double rhoU[], const double rhoW[], const double rhoTh[],
double rhoPrime[], double rhoUprime[], double rhoWprime[], double rhoThPrime[] ) {
    for( i = 0; i < M.N; i++ ) {
        P[i] = C.Po * pow( C.Rd * rhoTh[i] / C.Po, C.Cp / C.Cv );
    }
    //rho:
    singleVariableRHS( C, M, V, rho, rhoU, rhoW, rhoPrime );
    //rhoU:
    for( i = 0; i < M.N; i++ ) {
        F[i] = rhoU[i] * rhoU[i] / rho[i] + P[i];
        G[i] = rhoU[i] * rhoW[i] / rho[i];
    }
    singleVariableRHS( C, M, V, rhoU, F, G, rhoUprime );
    //rhoW:
    for( i = 0; i < M.N; i++ ) {
        F[i] = rhoW[i] * rhoW[i] / rho[i];
    }
    singleVariableRHS( C, M, V, rhoW, G, F, rhoWprime );
    for( k = 0; k < M.nLev; k++ ) { //adding missing terms ( - dP/dz - rho*g ):
        if( k == 0 ) { //lowest layer (extrapolate rho and enforce dP/dz = -rho*g):
            for( i = 0; i < M.n; i++ ) {
                tmpD = P[k*M.n+i] + C.g*M.dz * ( 3.*rho[k*M.n+i] - rho[(k+1)*M.n+i] ) / 2.;
                rhoWprime[k*M.n+i] = rhoWprime[k*M.n+i] - ( P[(k+1)*M.n+i] - tmpD ) / (2.*M.dz)
                - rho[k*M.n+i] * C.g;
            }
        }
        else if( k == M.nLev-1 ) { //highest layer (extrapolate rho and enforce dP/dz = -rho*g):
            for( i = 0; i < M.n; i++ ) {
                tmpD = P[k*M.n+i] - C.g*M.dz * ( 3.*rho[k*M.n+i] - rho[(k-1)*M.n+i] ) / 2.;
                rhoWprime[k*M.n+i] = rhoWprime[k*M.n+i] - ( tmpD - P[(k-1)*M.n+i] ) / (2.*M.dz)
                - rho[k*M.n+i] * C.g;
            }
        }
        else { //all other layers:
            for( i = 0; i < M.n; i++ ) {
                rhoWprime[k*M.n+i] = rhoWprime[k*M.n+i] - ( P[(k+1)*M.n+i] - P[(k-1)*M.n+i] ) / (2.*M.dz)
                - rho[k*M.n+i] * C.g;
            }
        }
    }
    //rhoTh:
    for( i = 0; i < M.N; i++ ) {
        F[i] = rhoU[i] * rhoTh[i] / rho[i];
        G[i] = rhoW[i] * rhoTh[i] / rho[i];
    }
    singleVariableRHS( C, M, V, rhoTh, F, G, rhoThPrime );
}

inline void TimeStepper::rk( const Constants& C, const FDmesh& M, Variables& V ) {
    if( rkStages == 2 ) {
        //2-stage, 2nd order RK.  The outputs are t and rho.  t increments and rho gets updated.
        //stage 1:
        odeFun( C, M, V, V.rho, V.rhoU, V.rhoW, V.rhoTh, s1_rho, s1_rhoU, s1_rhoW, s1_rhoTh );
        //stage 2:
        t = t + dt/2.;
        for( i = 0; i < M.N; i++ ) {
            tmp_rho[i]   = V.rho[i]   + dt/2. * s1_rho[i];
            tmp_rhoU[i]  = V.rhoU[i]  + dt/2. * s1_rhoU[i];
            tmp_rhoW[i]  = V.rhoW[i]  + dt/2. * s1_rhoW[i];
            tmp_rhoTh[i] = V.rhoTh[i] + dt/2. * s1_rhoTh[i];
        }
        odeFun( C, M, V, tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s1_rho, s1_rhoU, s1_rhoW, s1_rhoTh );
        //update t and get new value of rho:
        t = t + dt/2.;
        for( i = 0; i < M.N; i++ ) {
            V.rho[i]   = V.rho[i]   + dt * s1_rho[i];
            V.rhoU[i]  = V.rhoU[i]  + dt * s1_rhoU[i];
            V.rhoW[i]  = V.rhoW[i]  + dt * s1_rhoW[i];
            V.rhoTh[i] = V.rhoTh[i] + dt * s1_rhoTh[i];
            V.P[i] = C.Po * pow( C.Rd*V.rhoTh[i]/C.Po, C.Cp/C.Cv );
        }
    }
    else if( rkStages == 3 ) {
        //3-stage, 3rd order RK.  The outputs are t and rho.  t increments and rho gets updated.
        //stage 1:
        odeFun( C, M, V, V.rho, V.rhoU, V.rhoW, V.rhoTh, s1_rho, s1_rhoU, s1_rhoW, s1_rhoTh );
        //stage 2:
        t = t + dt/3.;
        for( i = 0; i < M.N; i++ ) {
            tmp_rho[i]   = V.rho[i]   + dt/3. * s1_rho[i];
            tmp_rhoU[i]  = V.rhoU[i]  + dt/3. * s1_rhoU[i];
            tmp_rhoW[i]  = V.rhoW[i]  + dt/3. * s1_rhoW[i];
            tmp_rhoTh[i] = V.rhoTh[i] + dt/3. * s1_rhoTh[i];
        }
        odeFun( C, M, V, tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s2_rho, s2_rhoU, s2_rhoW, s2_rhoTh );
        //stage 3:
        t = t + dt/3.;
        for( i = 0; i < M.N; i++ ) {
            tmp_rho[i]   = V.rho[i]   + 2*dt/3. * s2_rho[i];
            tmp_rhoU[i]  = V.rhoU[i]  + 2*dt/3. * s2_rhoU[i];
            tmp_rhoW[i]  = V.rhoW[i]  + 2*dt/3. * s2_rhoW[i];
            tmp_rhoTh[i] = V.rhoTh[i] + 2*dt/3. * s2_rhoTh[i];
        }
        odeFun( C, M, V, tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s2_rho, s2_rhoU, s2_rhoW, s2_rhoTh );
        //update t and get new value of rho:
        t = t + dt/3.;
        for( i = 0; i < M.N; i++ ) {
            V.rho[i]   = V.rho[i]   + dt/4. * ( s1_rho[i]   + 3*s2_rho[i] );
            V.rhoU[i]  = V.rhoU[i]  + dt/4. * ( s1_rhoU[i]  + 3*s2_rhoU[i] );
            V.rhoW[i]  = V.rhoW[i]  + dt/4. * ( s1_rhoW[i]  + 3*s2_rhoW[i] );
            V.rhoTh[i] = V.rhoTh[i] + dt/4. * ( s1_rhoTh[i] + 3*s2_rhoTh[i] );
            V.P[i] = C.Po * pow( C.Rd*V.rhoTh[i]/C.Po, C.Cp/C.Cv );
        } 
    }
    else if( rkStages == 4 ) {
        //4-stage, 4th order RK.  The outputs are t and rho.  t increments and rho gets updated.
        //stage 1:
        odeFun( C, M, V, V.rho, V.rhoU, V.rhoW, V.rhoTh, s1_rho, s1_rhoU, s1_rhoW, s1_rhoTh );
        //stage 2:
        t = t + dt/2.;
        for( i = 0; i < M.N; i++ ) {
            tmp_rho[i]   = V.rho[i]   + dt/2. * s1_rho[i];
            tmp_rhoU[i]  = V.rhoU[i]  + dt/2. * s1_rhoU[i];
            tmp_rhoW[i]  = V.rhoW[i]  + dt/2. * s1_rhoW[i];
            tmp_rhoTh[i] = V.rhoTh[i] + dt/2. * s1_rhoTh[i];
        }
        odeFun( C, M, V, tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s2_rho, s2_rhoU, s2_rhoW, s2_rhoTh );
        //stage 3:
        for( i = 0; i < M.N; i++ ) {
            tmp_rho[i]   = V.rho[i]   + dt/2. * s2_rho[i];
            tmp_rhoU[i]  = V.rhoU[i]  + dt/2. * s2_rhoU[i];
            tmp_rhoW[i]  = V.rhoW[i]  + dt/2. * s2_rhoW[i];
            tmp_rhoTh[i] = V.rhoTh[i] + dt/2. * s2_rhoTh[i];
        }
        odeFun( C, M, V, tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s3_rho, s3_rhoU, s3_rhoW, s3_rhoTh );
        //stage 4:
        t = t + dt/2;
        for( i = 0; i < M.N; i++ ) {
            tmp_rho[i]   = V.rho[i]   + dt * s3_rho[i];
            tmp_rhoU[i]  = V.rhoU[i]  + dt * s3_rhoU[i];
            tmp_rhoW[i]  = V.rhoW[i]  + dt * s3_rhoW[i];
            tmp_rhoTh[i] = V.rhoTh[i] + dt * s3_rhoTh[i];
        }
        odeFun( C, M, V, tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s4_rho, s4_rhoU, s4_rhoW, s4_rhoTh );
        //get new value:
        for( i = 0; i < M.N; i++ ) {
            V.rho[i]   = V.rho[i]   + dt/6. * ( s1_rho[i]   + 2*s2_rho[i]   + 2*s3_rho[i]   + s4_rho[i] );
            V.rhoU[i]  = V.rhoU[i]  + dt/6. * ( s1_rhoU[i]  + 2*s2_rhoU[i]  + 2*s3_rhoU[i]  + s4_rhoU[i] );
            V.rhoW[i]  = V.rhoW[i]  + dt/6. * ( s1_rhoW[i]  + 2*s2_rhoW[i]  + 2*s3_rhoW[i]  + s4_rhoW[i] );
            V.rhoTh[i] = V.rhoTh[i] + dt/6. * ( s1_rhoTh[i] + 2*s2_rhoTh[i] + 2*s3_rhoTh[i] + s4_rhoTh[i] );
            V.P[i] = C.Po * pow( C.Rd*V.rhoTh[i]/C.Po, C.Cp/C.Cv );
        }
    }
}

#endif