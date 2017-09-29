#include <iostream>
#include <cmath>

double r( double x, double y );
double phs( double r, double n );
double ga( double r, double ep );
double iq( double r, double ep );

int main()
{
	int m = 3;
	int n = 5;
	std::string rbfType("ga");
	double rbfParam = 1;
	double a = -1;
	double b = 1;
	double c = -1;
	double d = 1;
	
	double dx = ( b - a ) / ( n - 1 );
	double dy = ( d - c ) / ( m - 1 );
	double x[n];
	double y[m];
	double X[m*n];
	double Y[m*n];
	double A[m*n][m*n];
	
	std::cout << std::endl;
	
	//Get x vector:
	for( int i=0; i<n; i++ )
	{
		x[i] = a + i*dx;
		std::cout << "x[" << i << "] = " << x[i] << std::endl;
	}
	
	std::cout << std::endl;
	
	//Get y vector:
	for( int j=0; j<m; j++ )
	{
		y[j] = c + j*dy;
		std::cout << "y[" << j << "] = " << y[j] << std::endl;
	}
	
	std::cout << std::endl;
	
	//Get all x-coordinates and y coordinates (basically meshgrid):
	int k = 0;
	int rem = 0;
	X[0] = x[0];
	Y[0] = y[0];
	for( int i=1; i<m*n; i++ )
	{
		rem = i % m;
		if( rem == 0 )
		{
			k = k + 1;
		}
		X[i] = x[k];
		Y[i] = y[rem];
	}
	
	std::cout << std::endl;
	
	//Print x-coordinates
	for( int i=0; i<m*n; i++ )
	{
		std::cout << "X[" << i << "] = " << X[i] << std::endl;
	}
	
	std::cout << std::endl;
	
	//Print y-coordinates:
	for( int i=0; i<m*n; i++ )
	{
		std::cout << "Y[" << i << "] = " << Y[i] << std::endl;
	}
	
	std::cout << std::endl;
	
	//Get A-matrix:
	for( int i=0; i<m*n; i++ )
	{
		for( int j=0; j<m*n; j++ )
		{
			if( rbfType == "ga" )
			{
				A[i][j] = ga( r(X[i]-X[j],Y[i]-Y[j]), rbfParam );
			}
			else if( rbfType == "phs" )
			{
				A[i][j] = phs( r(X[i]-X[j],Y[i]-Y[j]), rbfParam );
			}
			else if( rbfType == "iq" )
			{
				A[i][j] = iq( r(X[i]-X[j],Y[i]-Y[j]), rbfParam );
			}
			std::cout << "A[" << i << "][" << j << "] = " << A[i][j] << std::endl;
		}
	}
}

double r( double x, double y )
{
	double r = pow( x*x + y*y, .5 );
	return r;
}

double phs( double r, double n )
{
	double z;
	if( ( std::abs(n-std::floor(n)) < 1e-14 ) && ( (int(n)%2) == 0 ) )
	{
		if( std::abs(r)<1e-14 )
		{
			z = 0;
		}
		else
		{
			z = pow(r,n) * std::log(r);
		}
	}
	else
	{
		z = pow(r,n);
	}
	return z;
}

double ga( double r, double ep )
{
	double z = std::exp( -(ep*r)*(ep*r) );
	return z;
}

double iq( double r, double ep )
{
	double z = 1 / ( 1 + (ep*r)*(ep*r) );
	return z;
}