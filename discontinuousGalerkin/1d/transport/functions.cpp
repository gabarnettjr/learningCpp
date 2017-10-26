//#include <iostream>

void getCardinalDerivatives4( const int& ne, const int& np, double x[], double dphi0dx[], double dphi1dx[], double dphi2dx[], double dphi3dx[] ) {
    //if( np != 4 ) {
    //    std::cerr << "Error:  This function should only be called if np=4." << std::endl;
    //}
    double x0;
    double x1;
    double x2;
    double x3;
    for( int i=0; i<ne; i++ ) {
        x0 = x[np*i];
        x1 = x[np*i+1];
        x2 = x[np*i+2];
        x3 = x[np*i+3];

        dphi0dx[np*i]   = ( (x0-x1)*(x0-x2) + (x0-x1)*(x0-x3) + (x0-x2)*(x0-x3) ) / ( (x0-x1)*(x0-x2)*(x0-x3) );
        dphi0dx[np*i+1] = (x1-x2)*(x1-x3) / ( (x0-x1)*(x0-x2)*(x0-x3) );
        dphi0dx[np*i+2] = (x2-x1)*(x2-x3) / ( (x0-x1)*(x0-x2)*(x0-x3) );
        dphi0dx[np*i+3] = (x3-x1)*(x3-x2) / ( (x0-x1)*(x0-x2)*(x0-x3) );

        dphi1dx[np*i]   = (x0-x2)*(x0-x3) / ( (x1-x0)*(x1-x2)*(x1-x3) );
        dphi1dx[np*i+1] = ( (x1-x0)*(x1-x2) + (x1-x0)*(x1-x3) + (x1-x2)*(x1-x3) ) / ( (x1-x0)*(x1-x2)*(x1-x3) );
        dphi1dx[np*i+2] = (x2-x0)*(x2-x3) / ( (x1-x0)*(x1-x2)*(x1-x3) );
        dphi1dx[np*i+3] = (x3-x0)*(x3-x2) / ( (x1-x0)*(x1-x2)*(x1-x3) );

        dphi2dx[np*i]   = (x0-x1)*(x0-x3) / ( (x2-x0)*(x2-x1)*(x2-x3) );
        dphi2dx[np*i+1] = (x1-x0)*(x1-x3) / ( (x2-x0)*(x2-x1)*(x2-x3) ); 
        dphi2dx[np*i+2] = ( (x2-x0)*(x2-x1) + (x2-x0)*(x2-x3) + (x2-x1)*(x2-x3) ) / ( (x2-x0)*(x2-x1)*(x2-x3) );
        dphi2dx[np*i+3] = (x3-x0)*(x3-x1) / ( (x2-x0)*(x2-x1)*(x2-x3) );

        dphi3dx[np*i]   = (x0-x1)*(x0-x2) / ( (x3-x0)*(x3-x1)*(x3-x2) );
        dphi3dx[np*i+1] = (x1-x0)*(x1-x2) / ( (x3-x0)*(x3-x1)*(x3-x2) );
        dphi3dx[np*i+2] = (x2-x0)*(x2-x1) / ( (x3-x0)*(x3-x1)*(x3-x2) );
        dphi3dx[np*i+3] = ( (x3-x0)*(x3-x1) + (x3-x0)*(x3-x2) + (x3-x1)*(x3-x2) ) / ( (x3-x0)*(x3-x1)*(x3-x2) );
    }
}

void odeFun( const int& ne, const int& np, const int& N, const double& u, double t, double w[], double rho[], double rhoPrime[], double dphi0dx[], double dphi1dx[], double dphi2dx[], double dphi3dx[] )
{
    for( int i=0; i<ne; i++ ) {
        if( i == 0 ) {
            rhoPrime[np*i] = ( u*rho[N-1] + u*rho[np*i] ) / 2. - u*( rho[np*i] - rho[N-1] );
        }
        else {
            rhoPrime[np*i] = ( u*rho[np*i-1] + u*rho[np*i] ) / 2. - u*( rho[np*i] - rho[np*i-1] );
        }
        rhoPrime[np*i+1] = 0;
        rhoPrime[np*i+2] = 0;
        if( i == ne-1 ) {
            rhoPrime[np*i+3] = -( ( u*rho[np*i+3] + u*rho[0] ) / 2. - u*( rho[0] - rho[np*i+3] ) );
        }
        else {
            rhoPrime[np*i+3] = -( ( u*rho[np*i+3] + u*rho[np*i+4] ) / 2. - u*( rho[np*i+4] - rho[np*i+3] ) );
        }
        if( np == 4 ) {
            for( int j=0; j<np; j++ ) {
                rhoPrime[np*i]   = rhoPrime[np*i]   + w[np*i+j]*u*rho[np*i+j]*dphi0dx[np*i+j];
                rhoPrime[np*i+1] = rhoPrime[np*i+1] + w[np*i+j]*u*rho[np*i+j]*dphi1dx[np*i+j];
                rhoPrime[np*i+2] = rhoPrime[np*i+2] + w[np*i+j]*u*rho[np*i+j]*dphi2dx[np*i+j];
                rhoPrime[np*i+3] = rhoPrime[np*i+3] + w[np*i+j]*u*rho[np*i+j]*dphi3dx[np*i+j];
            }
            rhoPrime[np*i]   = rhoPrime[np*i]   / w[np*i];
            rhoPrime[np*i+1] = rhoPrime[np*i+1] / w[np*i+1];
            rhoPrime[np*i+2] = rhoPrime[np*i+2] / w[np*i+2];
            rhoPrime[np*i+3] = rhoPrime[np*i+3] / w[np*i+3];
        }
    }
}

void rk( const int& ne, const int& np, const int& N, const double& u, double t, double w[], double rho[], double dphi0dx[], double dphi1dx[], double dphi2dx[], double dphi3dx[], const double& dt, double s1[], double s2[], double s3[], double s4[], double tmp[] ) {
    odeFun( ne, np, N, u, t,       w, rho, s1, dphi0dx, dphi1dx, dphi2dx, dphi3dx );
    for( int i=0; i<N; i++ ) {
        tmp[i] = rho[i] + dt/2.*s1[i];
    }
    odeFun( ne, np, N, u, t+dt/2., w, tmp, s2, dphi0dx, dphi1dx, dphi2dx, dphi3dx );
    for( int i=0; i<N; i++ ) {
        tmp[i] = rho[i] + dt/2.*s2[i];
    }
    odeFun( ne, np, N, u, t+dt/2., w, tmp, s3, dphi0dx, dphi1dx, dphi2dx, dphi3dx );
    for( int i=0; i<N; i++ ) {
        tmp[i] = rho[i] + dt*s3[i];
    }
    odeFun( ne, np, N, u, t+dt,    w, tmp, s4, dphi0dx, dphi1dx, dphi2dx, dphi3dx );
    for( int i=0; i<N; i++ ) {
        rho[i] = rho[i] + dt/6. * ( s1[i] + 2*s2[i] + 2*s3[i] + s4[i] );
    }
}

