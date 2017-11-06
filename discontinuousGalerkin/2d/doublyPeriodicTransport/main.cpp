#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "DGmesh.hpp"
#include "functions.hpp"

//Numerical solution for the 2d doubly periodic transport equation using discontinuous Galerkin
//laterally and second order finite differences vertically.

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
    const double c = -1.;                         //bottom endpoint
    const double d = 1.;                          //top endpoint
    const int np = 3;                             //number of polynomials per element
    const int ne = 20;                            //number of elements
    const int nLev = 40;                          //number of vertical levels
    const double dt = 1./100.;                    //time increment
    const int nTimesteps = 400;                   //number of timesteps
    const int rkStages = 3;                       //number of Runge-Kutta stages (2, 3, or 4)
    double t = 0.;                                //start time
    
    int i, j, k, ell;
    
    //Initialize mesh:
    DGmesh M( a, b, c, d, np, ne, nLev );
    
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
    
    //number of nodes per layer:
    const int n = M.getDFperLayer();
    
    //x-coordinates and quadrature weights, each given as a 1D array:
    double x[n];
    double weights[n];
    M.getCoordsAndQuadWeights( x, weights );
    
    //Create the cardinal derivatives:
    double dphi0dx[n];
    double dphi1dx[n];
    double dphi2dx[n];
    double dphi3dx[n];
    M.getCardinalDerivatives( dphi0dx, dphi1dx, dphi2dx, dphi3dx );
    
    //array of layer midpoints (z-coordinates):
    double z[nLev];
    M.getLayerMidpoints( z );
    
    //uniform layer thickness:
    const double dz = M.getLayerThickness();
    
    //total number of nodes (degrees of freedom):
    const int N = M.getDF();
    
    //initial density rho and velocity (u,w):
    double rho[N];
    double u[N];
    double w[N];
    for( ell=0; ell<nLev; ell++ ) {
        for( i=0; i<n; i++ ) {
            rho[ell*n+i] = rhoIC( x[i], z[ell] );
            u[ell*n+i] = uFunc( x[i], z[ell], t );
            w[ell*n+i] = wFunc( x[i], z[ell], t );
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
    
    //save layer thickness:
    outFile.open( "dz.txt" );
    outFile << dz;
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
    for( ell=0; ell<nTimesteps; ell++ ) {
        if( rkStages == 2 ) {
            rk2( i, j, k, ne, np, nLev, n, N, u, w, x, z, t, dz, alphaMax, weights, rho,
            dphi0dx, dphi1dx, dphi2dx, dphi3dx, dt, s1, s2, s3, s4, tmp, tmpD );
        }
        else if( rkStages == 3 ) {
            rk3( i, j, k, ne, np, nLev, n, N, u, w, x, z, t, dz, alphaMax, weights, rho,
            dphi0dx, dphi1dx, dphi2dx, dphi3dx, dt, s1, s2, s3, s4, tmp, tmpD );
        }
        else if( rkStages == 4 ) {
            rk4( i, j, k, ne, np, nLev, n, N, u, w, x, z, t, dz, alphaMax, weights, rho,
            dphi0dx, dphi1dx, dphi2dx, dphi3dx, dt, s1, s2, s3, s4, tmp, tmpD );
        }
        else {
            std::cerr << "Error:  rkStages should be 2, 3, or 4.";
            std::exit( EXIT_FAILURE );
        }
        std::stringstream s;
        s << "./snapshots/" << std::setfill('0') << std::setw(6) << ell+1 << ".txt";
        outFile.open( s.str() );
        for( i=0; i<N; i++ ) {
            outFile << rho[i] << " ";
        }
        outFile.close();
    }
    
    return 0;
}
