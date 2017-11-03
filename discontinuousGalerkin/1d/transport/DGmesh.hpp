#ifndef DGMESH_HPP
#define DGMESH_HPP

#include <cmath>

class DGmesh {
    
    public:
        
        DGmesh( double, double, int, int );
        int getN() { return N; }
        void getElementBoundaries( double[] );
        void getElementWidthsAndCenters( double[], double[] );
        void getGLL( double[], double[] );
        void getCoordsAndQuadWeights( double[], double[] );
        void getCardinalDerivatives( double[], double[], double[], double[] );
    
    private:
        
        double a;           //left endpoint
        double b;           //right endpoint
        int np;             //number of polynomials per element
        int ne;             //number of elements
        
        int N;              //total number of nodes (degrees of freedom)
        double* xGLL;       //GLL quadrature nodes on [-1,1]
        double* wGLL;       //GLL quadrature weights on [-1,1]
        double* xb;         //element boundaries
        double* dx;         //element widths
        double* xc;         //element centers
        double* x;          //all x-coordinates
        double* w;          //all quadrature weights
        double* dphi0dx;    //derivative of cardinal function phi0
        double* dphi1dx;    //derivative of cardinal function phi1
        double* dphi2dx;    //derivative of cardinal function phi2
        double* dphi3dx;    //derivative of cardinal function phi3
};

DGmesh::DGmesh( double A, double B, int NP, int NE ) {
    
    a = A;
    b = B;
    np = NP;
    ne = NE;
    
    N = np*ne;
    
    //element boundaries:
    xb = new double[ne+1];
    for( int i=0; i<ne+1; i++ ) {
        xb[i] = a + i*(b-a)/ne;
    }
    
    //element widths and centers:
    dx = new double[ne];
    xc = new double[ne];
    for( int i=0; i<ne; i++ ) {
        dx[i] = xb[i+1] - xb[i];
        xc[i] = ( xb[i] + xb[i+1] ) / 2.;
    }
    
    //GLL nodes and weights on [-1,1]:
    xGLL = new double[np];
    wGLL = new double[np];
    if( np == 2 ) {
        xGLL[0] = -1.;  xGLL[1] = 1.;
        wGLL[0] = 1.;  wGLL[1] = 1.;
    }
    else if( np == 3 ) {
        xGLL[0] = -1.;  xGLL[1] = 0.;  xGLL[2] = 1.;
        wGLL[0] = 1./3.;  wGLL[1] = 4./3.;  wGLL[2] = 1./3.;
    }
	else if( np == 4 ) {
        xGLL[0] = -1.;  xGLL[1] = -sqrt(5.)/5.;  xGLL[2] = sqrt(5.)/5;  xGLL[3] = 1.;
        wGLL[0] = 1./6.;  wGLL[1] = 5./6.;  wGLL[2] = 5./6.;  wGLL[3] = 1./6.;
    }
    
    //all nodes and quadrature weights:
    x = new double[N];
    w = new double[N];
    for( int i=0; i<ne; i++ ) {
        for( int j=0; j<np; j++ ) {
            x[np*i+j] = xc[i] + dx[i]/2. * xGLL[j];
            w[np*i+j] = dx[i]/2. * wGLL[j];
        }
    }
    
    //derivatives of cardinal functions:
    dphi0dx = new double[N];
    dphi1dx = new double[N];
    dphi2dx = new double[N];
    dphi3dx = new double[N];
    double x0, x1, x2, x3;
    for( int i=0; i<ne; i++ ) {
        if( np == 2 ) {
            x0 = x[np*i];
            x1 = x[np*i+1];
            
            dphi0dx[np*i]   = 1./(x0-x1);
            dphi0dx[np*i+1] = 1./(x0-x1);
            
            dphi1dx[np*i]   = 1./(x1-x0);
            dphi1dx[np*i+1] = 1./(x1-x0);
        }
        else if( np == 3) {
            x0 = x[np*i];
            x1 = x[np*i+1];
            x2 = x[np*i+2];
            
            dphi0dx[np*i]   = ( (x0-x1) + (x0-x2) ) / ( (x0-x1)*(x0-x2) );
            dphi0dx[np*i+1] = (x1-x2) / ( (x0-x1)*(x0-x2) );
            dphi0dx[np*i+2] = (x2-x1) / ( (x0-x1)*(x0-x2) );
            
            dphi1dx[np*i]   = (x0-x2) / ( (x1-x0)*(x1-x2) );
            dphi1dx[np*i+1] = ( (x1-x0) + (x1-x2) ) / ( (x1-x0)*(x1-x2) );
            dphi1dx[np*i+2] = (x2-x0) / ( (x1-x0)*(x1-x2) );
            
            dphi2dx[np*i]   = (x0-x1) / ( (x2-x0)*(x2-x1) );
            dphi2dx[np*i+1] = (x1-x0) / ( (x2-x0)*(x2-x1) );
            dphi2dx[np*i+2] = ( (x2-x0) + (x2-x1) ) / ( (x2-x0)*(x2-x1) );
        }
        else if( np == 4 ) {
            x0 = x[np*i];
            x1 = x[np*i+1];
            x2 = x[np*i+2];
            x3 = x[np*i+3];
            
            dphi0dx[np*i]   = ( (x0-x1)*(x0-x2) + (x0-x1)*(x0-x3) + (x0-x2)*(x0-x3) )
            / ( (x0-x1)*(x0-x2)*(x0-x3) );
            dphi0dx[np*i+1] = (x1-x2)*(x1-x3) / ( (x0-x1)*(x0-x2)*(x0-x3) );
            dphi0dx[np*i+2] = (x2-x1)*(x2-x3) / ( (x0-x1)*(x0-x2)*(x0-x3) );
            dphi0dx[np*i+3] = (x3-x1)*(x3-x2) / ( (x0-x1)*(x0-x2)*(x0-x3) );
            
            dphi1dx[np*i]   = (x0-x2)*(x0-x3) / ( (x1-x0)*(x1-x2)*(x1-x3) );
            dphi1dx[np*i+1] = ( (x1-x0)*(x1-x2) + (x1-x0)*(x1-x3) + (x1-x2)*(x1-x3) )
            / ( (x1-x0)*(x1-x2)*(x1-x3) );
            dphi1dx[np*i+2] = (x2-x0)*(x2-x3) / ( (x1-x0)*(x1-x2)*(x1-x3) );
            dphi1dx[np*i+3] = (x3-x0)*(x3-x2) / ( (x1-x0)*(x1-x2)*(x1-x3) );
            
            dphi2dx[np*i]   = (x0-x1)*(x0-x3) / ( (x2-x0)*(x2-x1)*(x2-x3) );
            dphi2dx[np*i+1] = (x1-x0)*(x1-x3) / ( (x2-x0)*(x2-x1)*(x2-x3) ); 
            dphi2dx[np*i+2] = ( (x2-x0)*(x2-x1) + (x2-x0)*(x2-x3) + (x2-x1)*(x2-x3) )
            / ( (x2-x0)*(x2-x1)*(x2-x3) );
            dphi2dx[np*i+3] = (x3-x0)*(x3-x1) / ( (x2-x0)*(x2-x1)*(x2-x3) );
            
            dphi3dx[np*i]   = (x0-x1)*(x0-x2) / ( (x3-x0)*(x3-x1)*(x3-x2) );
            dphi3dx[np*i+1] = (x1-x0)*(x1-x2) / ( (x3-x0)*(x3-x1)*(x3-x2) );
            dphi3dx[np*i+2] = (x2-x0)*(x2-x1) / ( (x3-x0)*(x3-x1)*(x3-x2) );
            dphi3dx[np*i+3] = ( (x3-x0)*(x3-x1) + (x3-x0)*(x3-x2) + (x3-x1)*(x3-x2) )
            / ( (x3-x0)*(x3-x1)*(x3-x2) );
        }
    }
}

//Accessors:

void DGmesh::getElementBoundaries( double XB[] ) {
    for( int i=0; i<ne+1; i++ ) {
        XB[i] = xb[i];
    }
}

void DGmesh::getElementWidthsAndCenters( double DX[], double XC[] ) {
    for( int i=0; i<ne; i++ ) {
        DX[i] = dx[i];
        XC[i] = xc[i];
    }
}

void DGmesh::getGLL( double Xgll[], double Wgll[] ) {
    for( int i=0; i<np; i++ ) {
        Xgll[i] = xGLL[i];
        Wgll[i] = wGLL[i];
    }
}

void DGmesh::getCoordsAndQuadWeights( double X[], double W[] ) {
    for( int i=0; i<N; i++ ) {
        X[i] = x[i];
        W[i] = w[i];
    }
}

void DGmesh::getCardinalDerivatives( double Dphi0dx[], double Dphi1dx[], double Dphi2dx[], double Dphi3dx[] ) {
    for( int i=0; i<N; i++ ) {
        Dphi0dx[i] = dphi0dx[i];
        Dphi1dx[i] = dphi1dx[i];
        Dphi2dx[i] = dphi2dx[i];
        Dphi3dx[i] = dphi3dx[i];
    }
}

#endif