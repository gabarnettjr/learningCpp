#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "functions.hpp"

//Numerical solution for the 1d periodic transport equation.

int main()
{
    //Parameters that the user chooses:
    const double a = -1.;                         //left endpoint
    const double b = 1.;                          //right endpoint
    const int np = 3;                             //number of polynomials per element
    const int ne = 40;                            //number of elements
    const double dt = 1./200.;                    //time increment
    const int nTimesteps = 400;                   //number of timesteps
    double t = 0.;                                //start time
    
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
    
    //initial rho and velocity u:
    double rho[N];
    double u[N];
    for( i=0; i<ne; i++ ) {
        for( j=0; j<np; j++ ) {
            rho[np*i+j] = rhoIC( x[np*i+j] );
            u[np*i+j] = uFunc( x[np*i+j], t );
        }
    }
    
    //Open filestream for saving things:
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
    
    //save vector of all x-coordinates:
    outFile.open( "x.txt" );
    for( i=0; i<N; i++ ) {
        outFile << x[i] << " ";
    }
    outFile.close();

    ////Save the velocity u:
    //outFile.open( "u.txt" );
    //for( i=0; i<N; i++ ) {
    //    outFile << u[i] << " ";
    //}
    //outFile.close();
    
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
    double alphaMax;
    double tmpD;
    for( k=0; k<nTimesteps; k++ ) {
        rk3( i, j, ne, np, N, u, x, t, alphaMax, w, rho, dphi0dx, dphi1dx, dphi2dx, dphi3dx,
        dt, s1, s2, s3, s4, tmp, tmpD );
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
