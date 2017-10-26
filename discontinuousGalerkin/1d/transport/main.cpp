#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

//This is a program to solve the 1d periodic transport equation.

//The initial condition for rho, which should be periodic:
double rhoIC( double& x ) {
    return exp( -10 * pow(x,2) );
    // return cos( M_PI * x );
}

void getCardinalDerivatives4( const int&, const int&, double[], double[], double[], double[], double[] );

void odeFun( int, int, const int&, const int&, const int&, const double&, double, double[], double[], double[], double[], double[], double[], double[] );

void rk( int, int, const int&, const int&, const int&, const double&, double,double[], double[], double[], double[], double[], double[], const double&, double[], double[], double[], double[], double[] );

int main()
{
    //Parameters that the user chooses:
    const double u = 1.;                          //velocity
    const double a = -1.;                         //left endpoint
    const double b = 1.;                          //right endpoint
    const int np = 4;                             //number of polynomials per element
    const int ne = 16;                            //number of elements
    double t = 0.;
    const double dt = 1./100.;
    const int nTimesteps = 200;
    
    int i, j, k;

    //Equally spaced element boundary points.  Alternatively,
    //these could come from some other place and they may
    //not be equally spaced.
    double xb[ne+1];
    for( i=0; i<ne+1; i++ ) {
        xb[i] = a + i*(b-a)/ne;
    }

    //The array of element widths (all the same on an equispaced grid):
    double dx[ne];
    for( i=0; i<ne; i++ ) {
        dx[i] = xb[i+1] - xb[i];
    }

    //Total number of nodes (element boundary nodes are repeated):
    int N = np*ne;

    //GLL nodes and weights on standard interval from -1 to 1:
    double xGLL[np];
    double wGLL[np];
    if( np == 2 ) {
        xGLL[0] = -1.;
        xGLL[1] = 1.;
        wGLL[0] = 1.;
        wGLL[1] = 1.;
    }
    else if( np == 3 ) {
        xGLL[0] = -1.;
        xGLL[1] = 0.;
        xGLL[2] = 1.;
        wGLL[0] = 1./3.;
        wGLL[1] = 4./3.;
        wGLL[2] = 1./3.;
    }
	else if( np == 4 ) {
        xGLL[0] = -1.;
        xGLL[1] = -sqrt(5.)/5.;
        xGLL[2] = sqrt(5.)/5;
        xGLL[3] = 1.;
        wGLL[0] = 1./6.;
        wGLL[1] = 5./6.;
        wGLL[2] = 5./6.;
        wGLL[3] = 1./6.;
    }

    //Center of mass of each element:
    double xc[ne];
    for( i=0; i<ne; i++ ) {
        xc[i] = ( xb[i] + xb[i+1] ) / 2.;
    }

    //All of the x-coordinates and quad weights in one long array.
    //Also all of initial values of rho:
    double x[N];
    double w[N];
    double rho[N];
    for( i=0; i<ne; i++ ) {
        for( j=0; j<np; j++ ) {
            x[np*i+j] = xc[i] + dx[i]/2. * xGLL[j];
            w[np*i+j] = dx[i]/2. * wGLL[j];
            rho[np*i+j] = rhoIC( x[np*i+j] );
        }
    }
    
    ////Print x and rho to make sure they are correct:
    //std::cout << std::endl;
    //for( i=0; i<N; i++ ) {
    //    std::cout << "x[" << i << "] = " << x[i] << std::endl;
    //}
    //std::cout << std::endl;
    //for( i=0; i<N; i++ ) {
    //    std::cout << "rho[" << i << "] = " << rho[i] << std::endl;
    //}
    //std::cout << std::endl;
    //for( i=0; i<N; i++ ) {
    //    std::cout << "w[" << i << "] = " << w[i] << std::endl;
    //}
    
    ////Check if integration is working:
    //std::cout << std::endl;
    //double I = 0;
    //for( i=0; i<N; i++ ) {
    //	I = I + w[i]*rho[i];
    //}
    //std::cout << "I = " << I << std::endl;

    //Create the cardinal derivatives:
    double dphi0dx[N];
    double dphi1dx[N];
    double dphi2dx[N];
    double dphi3dx[N];
    if( np == 4 ) {
        getCardinalDerivatives4( ne, np, x, dphi0dx, dphi1dx, dphi2dx, dphi3dx );
    }
    
    //Open filestream for saving things:
    std::ofstream outFile;
    outFile.open( "u.txt" );

    int pr = 16;

    //Save the constant velocity u:
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
    for( i=0; i<N-1; i++ ) {
        outFile << x[i] << " ";
    }
    outFile << x[N-1];
    outFile.close();

    //save vector of rho values at initial time:
    outFile.open( "./snapshots/00000.txt" );
    for( i=0; i<N-1; i++ ) {
        outFile << rho[i] << " ";
    }
    outFile << rho[N-1];
    outFile.close();

    std::stringstream s;
    s << "./snapshots/" <<  std::setfill('0') << std::setw(5) << 5 << ".txt";
    outFile.open( s.str() );
    outFile << rho[2];
    outFile.close();
    
    //Time stepping:
    double s1[N];
    double s2[N];
    double s3[N];
    double s4[N];
    double tmp[N];
    for( k=0; k<nTimesteps; k++ ) {
        rk( i, j, ne, np, N, u, t, w, rho, dphi0dx, dphi1dx, dphi2dx, dphi3dx, dt, s1, s2, s3, s4, tmp );
        t = t + dt;
        std::stringstream s;
        s << "./snapshots/" << std::setfill('0') << std::setw(5) << k+1 << ".txt";
        outFile.open( s.str() );
        for( i=0; i<N-1; i++ ) {
            outFile << rho[i] << " ";
        }
        outFile << rho[N-1];
        outFile.close();
    }
}
