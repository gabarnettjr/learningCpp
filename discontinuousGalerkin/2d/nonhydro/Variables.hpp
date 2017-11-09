#ifndef VARIABLES_HPP
#define VARIABLES_HPP

#include <cmath>

#include "Constants.hpp"
#include "DGmesh.hpp"

class Variables {
    
    public:
        
        Variables( const Constants&, const DGmesh& );
        
        double* rho;
        double* rhoU;
        double* rhoW;
        double* rhoTh;
        double* P;
        
        double* thetaBar;
        double* piBar;
        
    private:
        
        double* pi;
        double* th;
         
        void risingBubble( const Constants&, const DGmesh& );
};

Variables::Variables( const Constants& C, const DGmesh& M ) {
    rho      = new double[M.N];
    rhoU     = new double[M.N];
    rhoW     = new double[M.N];
    rhoTh    = new double[M.N];
    P        = new double[M.N];
    thetaBar = new double[M.N];
    piBar    = new double[M.N];
    pi       = new double[M.N];
    th       = new double[M.N];
    
    risingBubble( C, M );
}

inline void Variables::risingBubble( const Constants& C, const DGmesh& M ) {
    double R = 1500.;
    double xc = 5000.;
    double zc = 3000.;
    double k = 1./900.;
    double r;
    double thetaPrime0;
    for( int i = 0; i < M.nLev; i++ ) {
        for( int j = 0; j < M.n; j++ ) {
            thetaBar[i*M.n+j] = 300.;
            piBar[i*M.n+j] = 1. - C.g / C.Cp / thetaBar[i*M.n+j] * M.z[i];
            r = sqrt( pow(M.x[j]-xc,2.) + pow(M.z[i]-zc,2.) );
            thetaPrime0 = 2. * exp(-pow(k*r,2.));
            pi[i*M.n+j] = piBar[i*M.n+j];
            th[i*M.n+j] = thetaBar[i*M.n+j] + thetaPrime0;
        }
    }
    for( int i = 0; i < M.N; i++ ) {
        rho[i] = C.Po * pow( pi[i], C.Cv / C.Rd ) / C.Rd / th[i];
        rhoTh[i] = rho[i] * th[i];
        P[i] = C.Po * pow( C.Rd*rhoTh[i]/C.Po, C.Cp/C.Cv );
    }
    for( int i = 0; i < M.N; i++ ) {
        rhoU[i] = 0.;
        rhoW[i] = 0.;
    }
}

#endif