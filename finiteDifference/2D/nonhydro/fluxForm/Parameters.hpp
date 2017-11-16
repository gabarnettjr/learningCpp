#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <string>

class Parameters {
    
    public:
        
        Parameters( const std::string& );
        
        int stenX;              //stencil-size for centered lateral approximations (3 or 5)
        
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
    
    stenX = 5;
    
    if( s == "risingBubble" ) {
        a = 0.;
        b = 10000.;
        c = 0.;
        d = 10000.;
        n = 50;
        nLev = 100;
        rkStages = 3;
        t = 0.;
        dt = 1./4.;
        nTimesteps = 6000;
        saveDelta = 40;
    }
    else if( s == "inertiaGravityWaves" ) { //not working yet:
        a = 0.;
        b = 300000.;
        c = 0.;
        d = 10000.;
        n = 1200;
        nLev = 40;
        rkStages = 3;
        t = 0.;
        dt = 1./2.;
        nTimesteps = 6000;
        saveDelta = 20;
    }
    
    N = n * nLev;
    dx = (b-a)/n;
    dz = (d-c)/nLev;
}

#endif