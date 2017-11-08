#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "DG.hpp"

//Numerical solution for the 2d doubly periodic transport equation using
//discontinuous Galerkin laterally and second order finite differences vertically.

//Instructions:
// (*) mkdir snapshots
// (*) g++ functions.hpp DG.hpp main.cpp
// (*) ./a
// (*) python surfingScript.py (push any key to advance to next frame)

int main()
{   
    //Parameters that the user chooses:
    const double a = -1.;                         //left endpoint
    const double b = 1.;                          //right endpoint
    const double c = -1.;                         //bottom endpoint
    const double d = 1.;                          //top endpoint
    const int np = 4;                             //number of polynomials per element (2, 3, or 4)
    const int ne = 20;                            //number of elements per level (layer)
    const int nLev = 80;                          //number of levels (layers)
    const int rkStages = 4;                       //number of Runge-Kutta stages (2, 3, or 4)
    const int nTimesteps = 400;                   //number of timesteps
    double t = 0.;                                //start time
    const double dt = 1./100.;                    //time increment
    
    int i, j;
    
    //Initialize:
    DG M( a, b, c, d, np, ne, nLev, rkStages, nTimesteps, t, dt );
    
    //number of nodes per layer:
    const int n = M.getDFperLayer();
    
    //x-coordinates:
    double x[n];
    M.getXcoords( x );
    
    //array of layer midpoints (z-coordinates):
    double z[nLev];
    M.getLayerMidpoints( z );
    
    //total number of nodes (degrees of freedom):
    const int N = M.getDF();
    
    //initial density rho:
    double rho[N];
    for( j=0; j<nLev; j++ ) {
        for( i=0; i<n; i++ ) {
            rho[j*n+i] = rhoIC( x[i], z[j] );
        }
    }
    
    ///////////////////////////////////////////////////////////////////////
    
    //Open output file stream for saving things:
    std::ofstream outFile;
    
    //Precision for printing to text files:
    const int pr = 16;
    
    //Save number of polynomials per element:
    outFile.open( "np.txt" );
    outFile << np;
    outFile.close();
    
    //Save number of elements:
    outFile.open( "ne.txt" );
    outFile << ne;
    outFile.close();
    
    //Save start time:
    outFile.open( "t.txt" );
    outFile << std::scientific << std::setprecision(pr) << t;
    outFile.close();
    
    //Save time increment:
    outFile.open( "dt.txt" );
    outFile << dt;
    outFile.close();
    
    //Save number of time-steps:
    outFile.open( "nTimesteps.txt" );
    outFile << nTimesteps;
    outFile.close();
    
    //save array of x-coordinates:
    outFile.open( "x.txt" );
    for( i=0; i<n; i++ ) {
        outFile << x[i] << " ";
    }
    outFile.close();
    
    //save array of z-coordinates:
    outFile.open( "z.txt" );
    for( i=0; i<nLev; i++ ) {
        outFile << z[i] << " ";
    }
    outFile.close();
    
    //save array of rho values at initial time:
    outFile.open( "./snapshots/000000.txt" );
    for( i=0; i<N; i++ ) {
        outFile << rho[i] << " ";
    }
    outFile.close();
    
    ///////////////////////////////////////////////////////////////////////
    
    //Time stepping (and saving rho as it changes):
    for( j=0; j<nTimesteps; j++ ) {
        //advance rho and t with a single Runge-Kutta time step:
        M.rk( rho );
        //save new array rho:
        std::stringstream s;
        s << "./snapshots/" << std::setfill('0') << std::setw(6) << j+1 << ".txt";
        outFile.open( s.str() );
        for( i=0; i<N; i++ ) {
            outFile << rho[i] << " ";
        }
        outFile.close();
    }
    
    return 0;
}
