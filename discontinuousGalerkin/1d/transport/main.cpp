#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "functions.hpp"

//Numerical solution for the 1d periodic transport equation.

//Instructions:
//mkdir snapshots
//g++ functions.hpp main.cpp
//python plottingScript.py

//Function describing the initial condition for rho, which should be periodic:
double rhoIC( double& x ) {
    // return exp( -10 * pow(x,2) );
    return cos( M_PI * x );
    //if( ( x > -1./2. ) && ( x < 1./2. ) ) {
    //    return 1./2.;
    //}
    //else {
    //    return 0.;
    //}
}

int main()
{
    //Parameters that the user chooses:
    const double u = 1.;                          //constant velocity
    const double a = -1.;                         //left endpoint
    const double b = 1.;                          //right endpoint
    const int np = 3;                             //number of polynomials per element
    const int ne = 40;                            //number of elements
    double t = 0.;                                //start time
    const double dt = 1./200.;                    //time increment
    const int nTimesteps = 400;                   //number of timesteps
    
    int i, j, k;
    
    //Equally spaced element boundary points.  Alternatively, these
    //may come from some other place and might not be equally spaced:
    double xb[ne+1];
    getElementBoundaries( ne, a, b, xb );
    
    //The array of element widths (all the same on an equispaced grid),
    //and also the center of mass for each element:
    double dx[ne];
    double xc[ne];
    getElementWidthsAndCenters( ne, xb, dx, xc );
    
    //GLL nodes and weights on standard interval [-1,1]:
    double xGLL[np];
    double wGLL[np];
    i = getGLL( np, xGLL, wGLL );
    
    //x-coordinates and quadrature weights, each given as a 1D array:
    const int N = np*ne;
    double x[N];
    double w[N];
    getCoordsAndQuadWeights( ne, np, xc, dx, xGLL, wGLL, x, w );
    
    //Create the cardinal derivatives:
    double dphi0dx[N];
    double dphi1dx[N];
    double dphi2dx[N];
    double dphi3dx[N];
    i = getCardinalDerivatives( ne, np, x, dphi0dx, dphi1dx, dphi2dx, dphi3dx );
    
    //initial condition for rho:
    double rho[N];
    for( i=0; i<ne; i++ ) {
        for( j=0; j<np; j++ ) {
            rho[np*i+j] = rhoIC( x[np*i+j] );
        }
    }
    
    //Open filestream for saving things:
    std::ofstream outFile;
    
    //Precision for printing to text files:
    const int pr = 16;
    
    //Save the constant velocity u:
    outFile.open( "u.txt" );
    outFile << std::scientific << std::setprecision(pr) << u;
    outFile.close();
    
    //Save start time:
    outFile.open( "t.txt" );
    outFile << t;
    outFile.close();
    
    //Save time increment:
    outFile.open( "dt.txt" );
    outFile << dt;
    outFile.close();
    
    //Save number of time-steps:
    outFile.open( "nTimesteps.txt" );
    outFile << nTimesteps;
    outFile.close();
    
    //save vector of all x-coordinates:
    outFile.open( "x.txt" );
    for( i=0; i<N; i++ ) {
        outFile << x[i] << " ";
    }
    outFile.close();
    
    //save vector of rho values at initial time:
    outFile.open( "./snapshots/000000.txt" );
    for( i=0; i<N; i++ ) {
        outFile << rho[i] << " ";
    }
    outFile.close();
    
    //Time stepping:
    double s1[N];
    double s2[N];
    double s3[N];
    double s4[N];
    double tmp[N];
    for( k=0; k<nTimesteps; k++ ) {
        rk4( i, j, ne, np, N, u, t, w, rho, dphi0dx, dphi1dx, dphi2dx, dphi3dx, dt, s1, s2, s3, s4, tmp );
        t = t + dt;
        std::stringstream s;
        s << "./snapshots/" << std::setfill('0') << std::setw(6) << k+1 << ".txt";
        outFile.open( s.str() );
        for( i=0; i<N; i++ ) {
            outFile << rho[i] << " ";
        }
        outFile.close();
    }
    
    return 0;
}
