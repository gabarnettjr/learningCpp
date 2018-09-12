#ifndef TIMESTEPPER_HPP
#define TIMESTEPPER_HPP

#include <iostream>
#include <cmath>

#include "Parameters.hpp"
#include "Variables.hpp"

class TimeStepper {
    
    public:
        
        TimeStepper( const Parameters&, Variables& );
        
        void rk( const Parameters&, Variables& );
        
    private:
        
        double* p;          //temporary pressure for inside odefun
        double* s1_u;       //RK stage 1
        double* s1_w;
        double* s1_dpids;
        double* s1_th;
        double* s2_u;       //RK stage 2
        double* s2_w;
        double* s2_dpids;
        double* s2_th;
        double* s3_u;       //RK stage 3
        double* s3_w;
        double* s3_dpids;
        double* s3_th;
        double* s4_u;       //RK stage 4
        double* s4_w;
        double* s4_dpids;
        double* s4_th;
        double* tmp_u;      //used inside rk
        double* tmp_w;
        double* tmp_dpids;
        double* tmp_th;
        double tmpD;        //used inside odeFun
        double* F;
        double* G;
        int i, j, k;        //used in for-loops
        double* wx;         //first derivative weights
        double* whv;        //hyperviscosity weights
        double* u;          //temporary horizontal velocity array
        double* w;          //temporary vertical velocity array
        
        void singleVariableRHS( const Parameters&, const Variables&,
        const double[], const double[], const double[],
        double[] );
        
        void odeFun( const Parameters&, const Variables&,
        const double[], const double[], const double[], const double[],
        double[], double[], double[], double[] );
};

TimeStepper::TimeStepper( const Parameters& P, Variables& V ) {
    
    //arrays needed for RK time stepping:
    s1_u      = new double[P.N];
    s1_w      = new double[P.N];
    s1_dpids  = new double[P.N];
    s1_th     = new double[P.N];
    s2_u      = new double[P.N];
    s2_w      = new double[P.N];
    s2_dpids  = new double[P.N];
    s2_th     = new double[P.N];
    s3_u      = new double[P.N];
    s3_w      = new double[P.N];
    s3_dpids  = new double[P.N];
    s3_th     = new double[P.N];
    s4_u      = new double[P.N];
    s4_w      = new double[P.N];
    s4_dpids  = new double[P.N];
    s4_th     = new double[P.N];
    tmp_u     = new double[P.N];
    tmp_w     = new double[P.N];
    tmp_dpids = new double[P.N];
    tmp_th    = new double[P.N];
    
    //used inside singleVariableRHS:
    F = new double[P.N];
    G = new double[P.N];
    
    //arrays of finite difference weights:
    wx  = new double[P.stenX];
    whv = new double[P.stenX];
    if( P.stenX == 3 ) {
        wx[0] = -1./2.;
        wx[1] = 0.;
        wx[2] = 1./2.;
        whv[0] = 1.;
        whv[1] = -2.;
        whv[2] = 1.;
    }
    else if( P.stenX == 5 ) {
        wx[0] = 1./12.;
        wx[1] = -2./3.;
        wx[2] = 0.;
        wx[3] = 2./3.;
        wx[4] = -1./12.;
        whv[0] = 1.;
        whv[1] = -4.;
        whv[2] = 6.;
        whv[3] = -4.;
        whv[4] = 1.;
    }
    
    //temporary pressure and velocity variables:
    p = new double[P.N];
    u = new double[P.N];
    w = new double[P.N];
}

inline void TimeStepper::singleVariableRHS( const Parameters& P, const Variables& V,
const double rho[], const double rhoU[], const double rhoW[],
double rhoPrime[] ) {
    //Sets rhoPrime equal to -d(rhoU)/dx - d(rhoW)/dz, using FD for d/dx and d/dz.
    //This is used for each of the four prognostic variables: rho, rhoU, rhoW, rhoTh.
    //lateral operations using FD:
    if( P.stenX == 3 ) {
        for( k = 0; k < P.nLev; k++ ) {
            //i = 0:
            rhoPrime[k*P.n] = - ( wx[0]*rhoU[k*P.n+(P.n-1)] + wx[2]*rhoU[k*P.n+1] ) / P.dx
            + std::abs( u[k*P.n] ) * P.dx/2
            * ( whv[0]*rho[k*P.n+(P.n-1)] + whv[1]*rho[k*P.n] + whv[2]*rho[k*P.n+1] ) / pow(P.dx,2.);
            //mid-range i:
            for( i=1; i<P.n-1; i++ ) {
                rhoPrime[k*P.n+i] = - ( wx[0]*rhoU[k*P.n+(i-1)] + wx[2]*rhoU[k*P.n+(i+1)] ) / P.dx
                + std::abs( u[k*P.n+i] ) * P.dx/2
                * ( whv[0]*rho[k*P.n+(i-1)] + whv[1]*rho[k*P.n+i] + whv[2]*rho[k*P.n+(i+1)] ) / pow(P.dx,2.);
            }
            //i = P.n-1:
            rhoPrime[k*P.n+(P.n-1)] = - ( wx[0]*rhoU[k*P.n+(P.n-2)] + wx[2]*rhoU[k*P.n] ) / P.dx
            + std::abs( u[k*P.n+P.n-1] ) * P.dx/2
            * ( whv[0]*rho[k*P.n+(P.n-2)] + whv[1]*rho[k*P.n+(P.n-1)] + whv[2]*rho[k*P.n] ) / pow(P.dx,2.);
        }
    }
    else if( P.stenX == 5 ) {
        for( k = 0; k < P.nLev; k++ ) {
            //i = 0:
            rhoPrime[k*P.n] = - ( wx[0]*rhoU[k*P.n+(P.n-2)] + wx[1]*rhoU[k*P.n+(P.n-1)]
            + wx[3]*rhoU[k*P.n+1] + wx[4]*rhoU[k*P.n+2] ) / P.dx
            - std::abs( u[k*P.n] ) * pow(P.dx,3.)/12.
            * ( whv[0]*rho[k*P.n+(P.n-2)] + whv[1]*rho[k*P.n+(P.n-1)] + whv[2]*rho[k*P.n]
            + whv[3]*rho[k*P.n+1] + whv[4]*rho[k*P.n+1] ) / pow(P.dx,4.);
            //i = 1:
            rhoPrime[k*P.n+1] = - ( wx[0]*rhoU[k*P.n+(P.n-1)] + wx[1]*rhoU[k*P.n]
            + wx[3]*rhoU[k*P.n+2] + wx[4]*rhoU[k*P.n+3] ) / P.dx
            - std::abs( u[k*P.n+1] ) * pow(P.dx,3.)/12.
            * ( whv[0]*rho[k*P.n+(P.n-1)] + whv[1]*rho[k*P.n] + whv[2]*rho[k*P.n+1]
            + whv[3]*rho[k*P.n+2] + whv[4]*rho[k*P.n+3] ) / pow(P.dx,4.);
            //mid-range i:
            for( i = 2; i<P.n-2; i++ ) {
                rhoPrime[k*P.n+i] = - ( wx[0]*rhoU[k*P.n+(i-2)] + wx[1]*rhoU[k*P.n+(i-1)]
                + wx[3]*rhoU[k*P.n+(i+1)] + wx[4]*rhoU[k*P.n+(i+2)] ) / P.dx
                - std::abs( u[k*P.n+i] ) * pow(P.dx,3.)/12.
                * ( whv[0]*rho[k*P.n+(i-2)] + whv[1]*rho[k*P.n+(i-1)] + whv[2]*rho[k*P.n+i]
                + whv[3]*rho[k*P.n+(i+1)] + whv[4]*rho[k*P.n+(i+2)] ) / pow(P.dx,4.);
            }
            //i = P.n-2:
            rhoPrime[k*P.n+(P.n-2)] = - ( wx[0]*rhoU[k*P.n+(P.n-4)] + wx[1]*rhoU[k*P.n+(P.n-3)]
            + wx[3]*rhoU[k*P.n+(P.n-1)] + wx[4]*rhoU[k*P.n+(0)] ) / P.dx
            - std::abs( u[k*P.n+(P.n-2)] ) * pow(P.dx,3.)/12.
            * ( whv[0]*rho[k*P.n+(P.n-4)] + whv[1]*rho[k*P.n+(P.n-3)] + whv[2]*rho[k*P.n+(P.n-2)]
            + whv[3]*rho[k*P.n+(P.n-1)] + whv[4]*rho[k*P.n+(0)] ) / pow(P.dx,4.);
            //i = P.n-1:
            rhoPrime[k*P.n+(P.n-1)] = - ( wx[0]*rhoU[k*P.n+(P.n-3)] + wx[1]*rhoU[k*P.n+(P.n-2)]
            + wx[3]*rhoU[k*P.n+(0)] + wx[4]*rhoU[k*P.n+(1)] ) / P.dx
            - std::abs( u[k*P.n+(P.n-1)] ) * pow(P.dx,3.)/12.
            * ( whv[0]*rho[k*P.n+(P.n-3)] + whv[1]*rho[k*P.n+(P.n-2)] + whv[2]*rho[k*P.n+(P.n-1)]
            + whv[3]*rho[k*P.n+(0)] + whv[4]*rho[k*P.n+(0+1)] ) / pow(P.dx,4.);
        }
    }
    else {
        std::cerr << "Error: stenX should be 3 or 5.";
        std::exit( EXIT_FAILURE );
    }
    //vertical operations using FD:
    if( P.stenZ == 3 ) {
        for( i = 0; i < P.n; i++ ) {
            //k = 0:
            rhoPrime[i] = rhoPrime[i] - ( rhoW[P.n+i] + rhoW[i] ) / (2*P.dz)
            + std::abs( w[i] ) * P.dz/2.
            * ( rho[i] - 2.*rho[i] + rho[P.n+i] ) / pow(P.dz,2.);
            //mid-range k:
            for( k = 1; k < P.nLev-1; k++ ) {
                rhoPrime[k*P.n+i] = rhoPrime[k*P.n+i] - ( rhoW[(k+1)*P.n+i] - rhoW[(k-1)*P.n+i] ) / (2*P.dz)
                + std::abs( w[k*P.n+i] ) * P.dz/2.
                * ( rho[(k-1)*P.n+i] - 2.*rho[k*P.n+i] + rho[(k+1)*P.n+i] ) / pow(P.dz,2.);
            }
            //k = P.nLev-1:
            rhoPrime[(P.nLev-1)*P.n+i] = rhoPrime[(P.nLev-1)*P.n+i] - ( -rhoW[(P.nLev-1)*P.n+i] - rhoW[(P.nLev-2)*P.n+i] ) / (2*P.dz)
            + std::abs( w[(P.nLev-1)*P.n+i] ) * P.dz/2.
            * ( rho[(P.nLev-2)*P.n+i] - 2.*rho[(P.nLev-1)*P.n+i] + rho[(P.nLev-1)*P.n+i] ) / pow(P.dz,2.);
        }
    }
    else {
        std::cerr << "Error: stenZ should be 3.";
        std::exit( EXIT_FAILURE );
    }
}

inline void TimeStepper::odeFun( const Parameters& P, const Variables& V,
const double rho[], const double rhoU[], const double rhoW[], const double rhoTh[],
double rhoPrime[], double rhoUprime[], double rhoWprime[], double rhoThPrime[] ) {
    for( i = 0; i < P.N; i++ ) {
        p[i] = P.Po * pow( P.Rd*rhoTh[i]/P.Po, P.Cp/P.Cv );
        u[i] = rhoU[i] / rho[i];
        w[i] = rhoW[i] / rho[i];
    }
    //rho:
    singleVariableRHS( P, V, rho, rhoU, rhoW, rhoPrime );
    //rhoU:
    for( i = 0; i < P.N; i++ ) {
        F[i] = rhoU[i] * u[i] + p[i];
        G[i] = rhoU[i] * w[i];
    }
    singleVariableRHS( P, V, rhoU, F, G, rhoUprime );
    //rhoW:
    for( i = 0; i < P.N; i++ ) {
        F[i] = rhoW[i] * w[i];
    }
    singleVariableRHS( P, V, rhoW, G, F, rhoWprime );
    for( i = 0; i < P.n; i++ ) {
        //k = 0:
        tmpD = p[0*P.n+i] + P.g*P.dz * ( 3.*rho[0*P.n+i] - rho[(0+1)*P.n+i] ) / 2.;
        rhoWprime[0*P.n+i] = rhoWprime[0*P.n+i] - ( p[(0+1)*P.n+i] - tmpD ) / (2.*P.dz)
        - rho[0*P.n+i] * P.g;
        //mid-range k:
        for( k = 1; k < P.nLev-1; k++ ) {
            rhoWprime[k*P.n+i] = rhoWprime[k*P.n+i] - ( p[(k+1)*P.n+i] - p[(k-1)*P.n+i] ) / (2.*P.dz)
            - rho[k*P.n+i] * P.g;
        }
        //k = P.nLev-1:
        tmpD = p[(P.nLev-1)*P.n+i] - P.g*P.dz * ( 3.*rho[(P.nLev-1)*P.n+i] - rho[(P.nLev-2)*P.n+i] ) / 2.;
        rhoWprime[(P.nLev-1)*P.n+i] = rhoWprime[(P.nLev-1)*P.n+i] - ( tmpD - p[(P.nLev-2)*P.n+i] ) / (2.*P.dz)
        - rho[(P.nLev-1)*P.n+i] * P.g;
    }
    //rhoTh:
    for( i = 0; i < P.N; i++ ) {
        tmpD = rhoTh[i] / rho[i];
        F[i] = rhoU[i] * tmpD;
        G[i] = rhoW[i] * tmpD;
    }
    singleVariableRHS( P, V, rhoTh, F, G, rhoThPrime );
}

inline void TimeStepper::rk( const Parameters& P, Variables& V ) {
    if( P.rkStages == 2 ) {
        //2-stage, 2nd order RK.  The outputs are t and rho.  t increments and rho gets updated.
        //stage 1:
        odeFun( P, V, V.rho, V.rhoU, V.rhoW, V.rhoTh, s1_rho, s1_rhoU, s1_rhoW, s1_rhoTh );
        //stage 2:
        V.t = V.t + P.dt/2.;
        for( i = 0; i < P.N; i++ ) {
            tmp_rho[i]   = V.rho[i]   + P.dt/2. * s1_rho[i];
            tmp_rhoU[i]  = V.rhoU[i]  + P.dt/2. * s1_rhoU[i];
            tmp_rhoW[i]  = V.rhoW[i]  + P.dt/2. * s1_rhoW[i];
            tmp_rhoTh[i] = V.rhoTh[i] + P.dt/2. * s1_rhoTh[i];
        }
        odeFun( P, V, tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s1_rho, s1_rhoU, s1_rhoW, s1_rhoTh );
        //update t and get new value of rho:
        V.t = V.t + P.dt/2.;
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
        V.t = V.t + P.dt/3.;
        for( i = 0; i < P.N; i++ ) {
            tmp_rho[i]   = V.rho[i]   + P.dt/3. * s1_rho[i];
            tmp_rhoU[i]  = V.rhoU[i]  + P.dt/3. * s1_rhoU[i];
            tmp_rhoW[i]  = V.rhoW[i]  + P.dt/3. * s1_rhoW[i];
            tmp_rhoTh[i] = V.rhoTh[i] + P.dt/3. * s1_rhoTh[i];
        }
        odeFun( P, V, tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s2_rho, s2_rhoU, s2_rhoW, s2_rhoTh );
        //stage 3:
        V.t = V.t + P.dt/3.;
        for( i = 0; i < P.N; i++ ) {
            tmp_rho[i]   = V.rho[i]   + 2*P.dt/3. * s2_rho[i];
            tmp_rhoU[i]  = V.rhoU[i]  + 2*P.dt/3. * s2_rhoU[i];
            tmp_rhoW[i]  = V.rhoW[i]  + 2*P.dt/3. * s2_rhoW[i];
            tmp_rhoTh[i] = V.rhoTh[i] + 2*P.dt/3. * s2_rhoTh[i];
        }
        odeFun( P, V, tmp_rho, tmp_rhoU, tmp_rhoW, tmp_rhoTh, s2_rho, s2_rhoU, s2_rhoW, s2_rhoTh );
        //update t and get new value of rho:
        V.t = V.t + P.dt/3.;
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
        V.t = V.t + P.dt/2.;
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
        V.t = V.t + P.dt/2;
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