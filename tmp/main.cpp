#include <iostream>
#include <string>

void addTwoNumbers( const double& x, const double& y, double& z ) {
    z = x + y;
}

void addVecs( const int& n, double X[], double Y[], double Z[] ) {
    for( int i=0; i<n; i++ ) {
        Z[i] = X[i] + Y[i];
    }
}

int main()
{
    double x = 5.;
    double y =3.;
    double z;
    addTwoNumbers( x, y, z );
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
    std::cout << "z = " << z << std::endl << std::endl;
    addTwoNumbers( z, x, y );
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
    std::cout << "z = " << z << std::endl << std::endl;
    
    const int n = 5;
    double X[n];
    double Y[n];
    for( int i=0; i<n; i++ ) {
        X[i] = i;
        Y[i] = -i;
    }
    double Z[n];
    addVecs( n, X, Y, Z );
    for( int i=0; i<n; i++ ) {
        std::cout << "(X,Y,Z)(" << i << ") = ( " << X[i] << ", " << Y[i] << ", " << Z[i] << " )" << std::endl;
    }
}
