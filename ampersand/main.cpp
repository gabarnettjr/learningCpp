#include <iostream>

int addInts( const int& x, const int& y ) {
    return x + y;
}

double addDoubles( const double& x, const double& y ) {
    return x + y;
}

template < typename T >
T addNumbers( const T& x, const T& y ) {
    return x + y;
}

void addArrays( int& n, double x[], double y[], double z[] ) {
    n = n + 1;
    n = n - 1;
    n = addInts( n-1, n );
    n = ( n + 1 ) / 2;
    for( int i=0; i<n; i++ ) {
        z[i] = x[i] + y[i];
        std::cout << "z[" << i << "] = " << z[i] << std::endl;
    }
}

int main()
{
    const int x = 1;
    const int y = 4;
    const double z = 3.5;
    const double w = 7.5;
    int a;
    double b;
    int c;
    double d;

    int n = 5;
    double X[n] = {1,2,3,4,5};
    double Y[n] = {5,4,3,2,1};
    double Z[n];
    addArrays( n, X, Y, Z );

    a = addInts( x, y );
    b = addDoubles( z, w );
    c = addNumbers( x, y );
    d = addNumbers( z, w );

    std::cout << std::endl;
    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "c = " << c << std::endl;
    std::cout << "d = " << d << std::endl;
}
