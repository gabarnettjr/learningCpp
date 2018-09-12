#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <iostream>
#include <string>

class Parameters {
    
    public:
        
        Parameters( const std::string& );
        
        int stenX;              //stencil-size for centered lateral approximations (3 or 5)
        int stenZ;              //stencil-size for vertical approximations (only 3 works for now)
        
        double a;               //left endpoint
        double b;               //right endpoint
        double c;               //bottom endpoint
        double d;               //top endpoint
        int n;                  //degrees of freedom (nodes) per level
        int nLev;               //number of levels (layers)
        int rkStages;           //number of Runge-Kutta stages (2, 3, or 4)
        double t;               //start time
        double dt;              //time increment
        int nTimesteps;         //number of timesteps
        int saveDelta;          //number of time steps between saves  
        
        int N;                  //total number of degrees of freedom
        double dx;              //lateral spacing
        double dz;              //vertical spacing
        
        double Po = 100000.;    //reference pressure near bottom of domain
        double Rd = 287.04;     //gas constant for dry air
        double Cp = 1005.;      //specific heat at constant pressure
        double Cv = Cp - Rd;    //specific heat at constant volume
        double g = 9.81;        //gravitational acceleration
};

Parameters::Parameters( const std::string& s ) {
    t = 0;
    stenX = 5;
    stenZ = 3;
    rkStages = 3;
    if( s == "risingBubble" ) {
        a = 0.;
        b = 10000.;
        c = 0.;
        d = 10000.;
        n = int( (b-a) / 200. );
        nLev = int( (d-c) / 100. );
        dt = 1./4.;
        nTimesteps = 6000;
        saveDelta = 80;
    }
    else if( s == "inertiaGravityWaves" ) {
        a = 0.;
        b = 300000.;
        c = 0.;
        d = 10000.;
        n = int( (b-a) / 500. );
        nLev = int( (d-c) / 125. );
        dt = 1./3.;
        nTimesteps = 4500;
        saveDelta = 60;
    }
    else if( s == "densityCurrent" ) {
        a = -25600.;
        b = 25600.;
        c = 0.;
        d = 6400.;
        n = int( (b-a) / 100. );
        nLev = int( (d-c) / 100. );
        dt = 1./4.;
        nTimesteps = 3600;
        saveDelta = 40;
    }
    else if ( s == "movingDensityCurrent" ) {
        a = 0.;
        b = 36000;
        c = 0;
        d = 6400;
        n = int( (b-a) / 100. );
        nLev = int( (d-c) / 100. );
        dt = 1./6.;
        nTimesteps = 5400;
        saveDelta = 60;
    }
    else {
        std::cerr << "Error: Invalid test case string s.";
        std::exit( EXIT_FAILURE );
    }
    N = n * nLev;
    dx = (b-a)/n;
    dz = (d-c)/nLev;
}

#endif