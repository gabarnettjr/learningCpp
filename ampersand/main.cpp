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
