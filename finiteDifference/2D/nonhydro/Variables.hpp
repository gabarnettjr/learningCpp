#ifndef VARIABLES_HPP
#define VARIABLES_HPP

#include <iostream>
#include <cmath>
#include <string>

#include "Parameters.hpp"
#include "FDmesh.hpp"

class Variables {
    
    public:
        
        Variables( const std::string&, const Parameters&, const FDmesh& );
        
        double t;           //time
        
        double* rho;        //density
        double* rhoU;       //density * ( horizontal velocity )
        double* rhoW;       //density * ( vertical velocity )
        double* rhoTh;      //density * ( potential temperature )
        double* p;          //pressure
        
        double* rhoBar;     //background density
        double* thetaBar;   //background potential temperature
        double* piBar;      //background Exner pressure
        double* pBar;       //background pressure
        
    private:
        
        double* pi;
        double* th;
         
        void risingBubble( const Parameters&, const FDmesh& );
        void densityCurrent( const Parameters&, const FDmesh& );
        void movingDensityCurrent( const Parameters&, const FDmesh& );
        void inertiaGravityWaves( const Parameters&, const FDmesh& );
};

Variables::Variables( const std::string& s, const Parameters& P, const FDmesh& M ) {
    t = P.t;
    
    rho      = new double[P.N];
    rhoU     = new double[P.N];
    rhoW     = new double[P.N];
    rhoTh    = new double[P.N];
    p        = new double[P.N];
    
    rhoBar   = new double[P.N];
    thetaBar = new double[P.N];
    piBar    = new double[P.N];
    pBar     = new double[P.N];
    
    pi       = new double[P.N];
    th       = new double[P.N];
    
    if( s == "risingBubble" ) {
        risingBubble( P, M );
    }
    else if( s == "inertiaGravityWaves" ) {
        inertiaGravityWaves( P, M );
    }
    else if( s == "densityCurrent" ) {
        densityCurrent( P, M );
    }
    else if( s == "movingDensityCurrent" ) {
        movingDensityCurrent( P, M );
    }
    else {
        std::cerr << "Error: Invalid test case string s.";
        std::exit( EXIT_FAILURE );
    }
}

inline void Variables::risingBubble( const Parameters& P, const FDmesh& M ) {
    double R = 1500.;
    double xc = 5000.;
    double zc = 3000.;
    double k = 1./900.;
    double r;
    double thetaPrime0;
    for( int i = 0; i < P.nLev; i++ ) {
        for( int j = 0; j < P.n; j++ ) {
            thetaBar[i*P.n+j] = 300.;
            piBar[i*P.n+j] = 1. - P.g / P.Cp / thetaBar[i*P.n+j] * M.z[i];
            r = sqrt( pow(M.x[j]-xc,2.) + pow(M.z[i]-zc,2.) );
            thetaPrime0 = 2. * exp(-pow(k*r,2.));
            th[i*P.n+j] = thetaBar[i*P.n+j] + thetaPrime0;
        }
    }
    for( int i = 0; i < P.N; i++ ) {
        pi[i] = piBar[i];
        rhoBar[i] = P.Po * pow( piBar[i], P.Cv/P.Rd ) / P.Rd / thetaBar[i];
        pBar[i] = P.Po * pow( P.Rd * rhoBar[i] * thetaBar[i] / P.Po, P.Cp/P.Cv );
        rho[i] = P.Po * pow( pi[i], P.Cv / P.Rd ) / P.Rd / th[i];
        rhoTh[i] = rho[i] * th[i];
        p[i] = P.Po * pow( P.Rd*rhoTh[i]/P.Po, P.Cp/P.Cv );
        rhoU[i] = 0.;
        rhoW[i] = 0.;
    }
}

inline void Variables::inertiaGravityWaves( const Parameters& P, const FDmesh& M ) {
    double N = .01;
    double theta0 = 300.;
    double thetaC = .01;
    double hC = 10000.;
    double aC = 5000.;
    double xC = 100000.;
    double thetaPrime0;
    for( int i = 0; i < P.nLev; i++ ) {
        for( int j = 0; j < P.n; j++ ) {
            thetaBar[i*P.n+j] = theta0 * exp( (pow(N,2.)/P.g) * M.z[i] );
            piBar[i*P.n+j] = 1. + pow(P.g,2.) / P.Cp / theta0 / pow(N,2) * ( exp(-pow(N,2.)/P.g*M.z[i]) - 1. );
            // piBar[i*P.n+j] = (1.-pow(P.g,2.)/pow(N,2.)/P.Cp/theta0) + pow(P.g,2.)/pow(N,2.)/P.Cp/theta0*exp(-(pow(N,2.)/P.g)*M.z[i]);
            thetaPrime0 = thetaC * sin( M_PI*M.z[i]/hC ) / ( 1. + pow((M.x[j]-xC)/aC,2.) );
            th[i*P.n+j] = thetaBar[i*P.n+j] + thetaPrime0;
        }
    }
    for( int i = 0; i < P.N; i++ ) {
        pi[i] = piBar[i];
        rhoBar[i] = P.Po * pow( piBar[i], P.Cv/P.Rd ) / P.Rd / thetaBar[i];
        pBar[i] = P.Po * pow( P.Rd * rhoBar[i] * thetaBar[i] / P.Po, P.Cp/P.Cv );
        rho[i] = P.Po * pow( pi[i], P.Cv / P.Rd ) / P.Rd / th[i];
        rhoTh[i] = rho[i] * th[i];
        p[i] = P.Po * pow( P.Rd*rhoTh[i]/P.Po, P.Cp/P.Cv );
        rhoU[i] = rho[i] * 20.;
        rhoW[i] = 0.;
    }
}

inline void Variables::densityCurrent( const Parameters& P, const FDmesh& M ) {
    double xc = 0.;
    double zc = 3000.;
    double xr = 4000.;
    double zr = 2000.;
    double rTilde;
    double Tprime0;
    double thetaPrime0;
    for( int i = 0; i < P.nLev; i++ ) {
        for( int j = 0; j < P.n; j++ ) {
            thetaBar[i*P.n+j] = 300.;
            piBar[i*P.n+j] = 1. - P.g / P.Cp / thetaBar[i*P.n+j] * M.z[i];
            rTilde = sqrt( pow((M.x[j]-xc)/xr,2.) + pow((M.z[i]-zc)/zr,2.) );
            if( rTilde <= 1. ) {
                Tprime0 = -15./2. * ( 1. + cos(M_PI*rTilde) );
            }
            else {
                Tprime0 = 0.;
            }
            thetaPrime0 = Tprime0 / piBar[i*P.n+j];
            th[i*P.n+j] = thetaBar[i*P.n+j] + thetaPrime0;
        }
    }
    for( int i = 0; i < P.N; i++ ) {
        pi[i] = piBar[i];
        rhoBar[i] = P.Po * pow( piBar[i], P.Cv/P.Rd ) / P.Rd / thetaBar[i];
        pBar[i] = P.Po * pow( P.Rd * rhoBar[i] * thetaBar[i] / P.Po, P.Cp/P.Cv );
        rho[i] = P.Po * pow( pi[i], P.Cv / P.Rd ) / P.Rd / th[i];
        rhoTh[i] = rho[i] * th[i];
        p[i] = P.Po * pow( P.Rd*rhoTh[i]/P.Po, P.Cp/P.Cv );
        rhoU[i] = 0.;
        rhoW[i] = 0.;
    }
}

inline void Variables::movingDensityCurrent( const Parameters& P, const FDmesh& M ) {
    double xc = 18000.;
    double zc = 3000.;
    double xr = 4000.;
    double zr = 2000.;
    double rTilde;
    double Tprime0;
    double thetaPrime0;
    for( int i = 0; i < P.nLev; i++ ) {
        for( int j = 0; j < P.n; j++ ) {
            thetaBar[i*P.n+j] = 300.;
            piBar[i*P.n+j] = 1. - P.g / P.Cp / thetaBar[i*P.n+j] * M.z[i];
            rTilde = sqrt( pow((M.x[j]-xc)/xr,2.) + pow((M.z[i]-zc)/zr,2.) );
            if( rTilde <= 1. ) {
                Tprime0 = -15./2. * ( 1. + cos(M_PI*rTilde) );
            }
            else {
                Tprime0 = 0.;
            }
            thetaPrime0 = Tprime0 / piBar[i*P.n+j];
            th[i*P.n+j] = thetaBar[i*P.n+j] + thetaPrime0;
        }
    }
    for( int i = 0; i < P.N; i++ ) {
        pi[i] = piBar[i];
        rhoBar[i] = P.Po * pow( piBar[i], P.Cv/P.Rd ) / P.Rd / thetaBar[i];
        pBar[i] = P.Po * pow( P.Rd * rhoBar[i] * thetaBar[i] / P.Po, P.Cp/P.Cv );
        rho[i] = P.Po * pow( pi[i], P.Cv / P.Rd ) / P.Rd / th[i];
        rhoTh[i] = rho[i] * th[i];
        p[i] = P.Po * pow( P.Rd*rhoTh[i]/P.Po, P.Cp/P.Cv );
        rhoU[i] = rho[i] * 20.;
        rhoW[i] = 0.;
    }
}

#endif