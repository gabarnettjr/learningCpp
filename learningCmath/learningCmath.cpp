#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>

void meshgridCols( double x[], double y[], const int& m, const int& n, double XY[][2] );
void r( double x, double y, double& rij );
double phs( double& r, double& n );
double ga( double& r, double& ep );
double iq( double& r, double& ep );
double mq( double& r, double& ep );

int main( int argc, char* argv[] )
{
	std::string rbfType( argv[1] );
	double rbfParam;  std::stringstream(argv[2]) >> rbfParam;
	int m;  std::stringstream(argv[3]) >> m;
	int n;  std::stringstream(argv[4]) >> n;

	const double a = -1;
	const double b = 1;
	const double c = -1;
	const double d = 1;
	const double dx = ( b - a ) / ( n - 1 );
	const double dy = ( d - c ) / ( m - 1 );
	
	double x[n];
	double y[m];
	double XY[m*n][2];
	double A[m*n][m*n];
	double rij;
	int i;
	int j;

	// std::cout << std::endl;
	
	//Get x vector:
	for( i=0; i<n; i++ )
	{
		x[i] = a + i*dx;
		// std::cout << "x[" << i << "] = " << x[i] << std::endl;
	}
	
	// std::cout << std::endl;
	
	//Get y vector:
	for( j=0; j<m; j++ )
	{
		y[j] = c + j*dy;
		// std::cout << "y[" << j << "] = " << y[j] << std::endl;
	}
	
	// std::cout << std::endl;
	
	//Get all x-coordinates and y coordinates (basically meshgrid):
	meshgridCols( x, y, m, n, XY );
	
	// std::cout << std::endl;
	
	//Print x-coordinates
	// for( int i=0; i<m*n; i++ )
	// {
		// std::cout << "X[" << i << "] = " << X[i] << std::endl;
	// }
	
	// std::cout << std::endl;
	
	//Print y-coordinates:
	// for( int i=0; i<m*n; i++ )
	// {
		// std::cout << "Y[" << i << "] = " << Y[i] << std::endl;
	// }
	
	// std::cout << std::endl;
	
	//Get A-matrix:
	for( i=0; i<m*n; i++ )
	{
		for( j=0; j<m*n; j++ )
		{
			r( XY[i][0] - XY[j][0], XY[i][1] - XY[j][1], rij );
			if( rbfType == "ga" )
			{
				A[i][j] = ga( rij, rbfParam );
			}
			else if( rbfType == "phs" )
			{
				A[i][j] = phs( rij, rbfParam );
			}
			else if( rbfType == "iq" )
			{
				A[i][j] = iq( rij, rbfParam );
			}
			else if( rbfType == "mq" )
			{
				A[i][j] = mq( rij, rbfParam );
			}
			else
			{
				std::cerr << "Failure.  Pick a valid RBF acronym." << std::endl;
				return(EXIT_FAILURE);
			}
			std::cout << "A[" << i << "][" << j << "] = " << A[i][j] << std::endl;
		}
	}
}

void meshgridCols( double x[], double y[], const int& m, const int& n, double XY[][2] )
{
	int k = 0;
	int rem = 0;
	int i;
	XY[0][0] = x[0];
	XY[0][1] = y[0];
	for( i=1; i<m*n; i++ )
	{
		rem = i % m;
		if( rem == 0 )
		{
			k = k + 1;
		}
		XY[i][0] = x[k];
		XY[i][1] = y[rem];
	}
}

void r( double x, double y, double& rij )
{
	rij = pow( x*x + y*y, .5 );
}

//void( double rbfType, double XY[][2], double A[][int mn] )
//{
	
//}

double phs( double& r, double& n )
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

double ga( double& r, double& ep )
{
	double z = std::exp( -(ep*r)*(ep*r) );
	return z;
}

double iq( double& r, double& ep )
{
	double z = 1 / ( 1 + (ep*r)*(ep*r) );
	return z;
}

double mq( double& r, double& ep )
{
	double z = std::sqrt( 1 + (ep*r)*(ep*r) );
	return z;
}
