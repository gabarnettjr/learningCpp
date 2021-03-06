#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "DGmesh.hpp"
#include "functions.hpp"

//Numerical solution for the 1d periodic transport equation using discontinuous Galerkin.

//Instructions:
// (*) mkdir snapshots
// (*) g++ DGmesh.hpp functions.hpp main.cpp
// (*) ./a
// (*) python plottingScript.py (push any key to advance to next frame)

int main()
{   
    //Parameters that the user chooses:
    const double a = -1.;                         //left endpoint
    const double b = 1.;                          //right endpoint
    const int np = 3;                             //number of polynomials per element
    const int ne = 40;                            //number of elements
    const double dt = 1./200.;                    //time increment
    const int nTimesteps = 800;                   //number of timesteps
    const int rkStages = 3;                       //number of Runge-Kutta stages (2, 3, or 4)
    double t = 0.;                                //start time
    
    int i, j, k;
    
    //Initialize mesh:
    DGmesh M( a, b, np, ne );
    
    //Equally spaced element boundary points.  Alternatively, these
    //may come from some other place and might not be equally spaced:
    double xb[ne+1];
    M.getElementBoundaries( xb );
    
    //The array of element widths (all the same on an equispaced grid),
    //and also the center of mass for each element:
    double dx[ne];
    double xc[ne];
    M.getElementWidthsAndCenters( dx, xc );
    
    //GLL nodes and weights on standard interval [-1,1]:
    double xGLL[np];
    double wGLL[np];
    M.getGLL( xGLL, wGLL );
    
    //total number of nodes (degrees of freedom)
    const int N = M.getN();
    
    //x-coordinates and quadrature weights, each given as a 1D array:
    double x[N];
    double w[N];
    M.getCoordsAndQuadWeights( x, w );
    
    //Create the cardinal derivatives:
    double dphi0dx[N];
    double dphi1dx[N];
    double dphi2dx[N];
    double dphi3dx[N];
    M.getCardinalDerivatives( dphi0dx, dphi1dx, dphi2dx, dphi3dx );
    
    //initial rho and velocity u:
    double rho[N];
    double u[N];
    for( i=0; i<ne; i++ ) {
        for( j=0; j<np; j++ ) {
            rho[np*i+j] = rhoIC( x[np*i+j] );
            u[np*i+j] = uFunc( x[np*i+j], t );
        }
    }
    
    //Open output file stream for saving things:
    std::ofstream outFile;
    
    //Precision for printing to text files:
    const int pr = 16;
    
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
    
    //save array of all x-coordinates:
    outFile.open( "x.txt" );
    for( i=0; i<N; i++ ) {
        outFile << x[i] << " ";
    }
    outFile.close();
    
    //save array of rho values at initial time:
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
    double alphaMax;
    double tmpD;
    for( k=0; k<nTimesteps; k++ ) {
        if( rkStages == 2 ) {
            rk2( i, j, ne, np, N, u, x, t, alphaMax, w, rho, dphi0dx, dphi1dx, dphi2dx, dphi3dx,
            dt, s1, s2, s3, s4, tmp, tmpD );
        }
        else if( rkStages == 3 ) {
            rk3( i, j, ne, np, N, u, x, t, alphaMax, w, rho, dphi0dx, dphi1dx, dphi2dx, dphi3dx,
            dt, s1, s2, s3, s4, tmp, tmpD );
        }
        else if( rkStages == 4 ) {
            rk4( i, j, ne, np, N, u, x, t, alphaMax, w, rho, dphi0dx, dphi1dx, dphi2dx, dphi3dx,
            dt, s1, s2, s3, s4, tmp, tmpD );
        }
        else {
            std::cerr << "Error:  rkStages should be 2, 3, or 4.";
            std::exit( EXIT_FAILURE );
        }
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
