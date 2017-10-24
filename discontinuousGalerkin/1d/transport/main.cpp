#include <iostream>
#include <cmath>

//This is a program to solve the 1d periodic transport equation.

//The initial condition for rho, which should be periodic:
double rhoIC( double x ) {
    return exp( -10 * pow(x,2) );
    // return cos( M_PI * x );
}

//double weights( int n, double z[], double x[], int m, double c[] ) {
//    double c1 = 1.;
//    double c4;
//    for( int i=0; i<n; i++ ) {
//        c4[i] = x[0] - z[i];
//    }
//    for( int i=1; i<n; i++ ) {
//        
//    }
//}
//
//void odeFun( int ne, int np, double u, double t, double rho[], double rhoPrime[] )
//{
//    for( int i=0; i<ne; i++ ) {
//        rhoPrime[np*i] = (  ) / 2. - u*(  ) + 
//        for( int j=1; j<np-1; j++ ) {
//            rhoPrime[np*i+j] = ;
//        }
//        rhoPrime[np*i+np-1] = ;
//    }
//}
//
//void rk( int ne, int np, double u, double t, double rho, double s1[], double s2[], double s3[], double s4[] ) {
//    odeFun( ne, np, u, t, rho, s1 );
//    odeFun( ne, np, u, t+dt/2, rho+dt/2*s1, s2 );
//    odeFun( ne, np, u, t+dt/2, rho+dt/2*s2, s3 );
//    odeFun( ne, np, u, t+dt, rho+dt*s3, s4 );
//    rho = rho + dt/6 * ( s1 + 2*s2 + 2*s3 + s4 );
//}

int main()
{
    //Parameters that the user chooses:
    double u = 1.;                          //velocity
    double a = -1.;                         //left endpoint
    double b = 1.;                          //right endpoint
    int np = 3;                             //number of polynomials per element
    int ne = 32;                             //number of elements
    
    int i, j, k;

    //Equally spaced element boundary points.  Alternatively,
    //these could come from some other place and they may
    //not be equally spaced.
    double x[ne+1];
    for( i=0; i<ne+1; i++ ) {
        x[i] = a + i*(b-a)/ne;
    }

    //The array of element widths (all the same on an equispaced grid):
    double dx[ne];
    for( i=0; i<ne; i++ ) {
        dx[i] = x[i+1] - x[i];
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
        xc[i] = ( x[i] + x[i+1] ) / 2.;
    }

    //All of the x-coordinates and quad weights in one long array.
    //Also all of initial values of rho:
    double X[N];
    double W[N];
    double rho[N];
    for( i=0; i<ne; i++ ) {
        for( j=0; j<np; j++ ) {
            X[np*i+j] = xc[i] + dx[i]/2. * xGLL[j];
            W[np*i+j] = dx[i]/2. * wGLL[j];
            rho[np*i+j] = rhoIC( X[np*i+j] );
        }
    }
    
    ////Print X and rho to make sure they are correct:
    //std::cout << std::endl;
    //for( i=0; i<N; i++ ) {
    //    std::cout << "X[" << i << "] = " << X[i] << std::endl;
    //}
    //std::cout << std::endl;
    //for( i=0; i<N; i++ ) {
    //    std::cout << "rho[" << i << "] = " << rho[i] << std::endl;
    //}
    //std::cout << std::endl;
    //for( i=0; i<N; i++ ) {
    //    std::cout << "W[" << i << "] = " << W[i] << std::endl;
    //}
    
    //Check if integration is working:
    std::cout << std::endl;
    double I = 0;
    for( i=0; i<N; i++ ) {
    	I = I + W[i]*rho[i];
    }
    std::cout << "I = " << I << std::endl;
    
    ////Time stepping:
    //double s1[N];
    //double s2[N];
    //double s3[N];
    //double s4[N];
    //for( k=0; k<nTimesteps; k++ ) {
    //    rk( ne, np, u, t, rho, s1, s2, s3, s4 );
    //    t = t + dt;
    //}
}
