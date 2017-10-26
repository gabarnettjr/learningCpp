#include <iostream>
//inserting this comment to test in git.
//inserting one more comment for one more test.

void addNums( const double&, const double&, double& );

void addVecs( const int&, double[], double[], double[] );

int main()
{
    double x = 5.;
    double y =3.;
    double z;
    addNums( x, y, z );
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
    std::cout << "z = " << z << std::endl << std::endl;
    addNums( z, x, y );
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
