#include <iostream>
#include <cmath>
#include <vector>

int main( int argc, char * argv[] )
{
	{
    	double x = exp(-2);
    	std::cout << "x = " << x << std::endl;
	}
	{
		int x = 4;
		std::cout << "x = " << x << std::endl;
	}
	{
		int n = 10;
		std::vector<double> x;
		for( int i=0; i<n; i++ )
		{
			x.push_back( i*i );
		}
    	for( int i=0; i<n; i++ )
		{
    	    std::cout << "x[" <<  i << "] = " << x[i] << std::endl;
    	}
		{
			std::vector<double> y;
			std::vector<double> z;
			for( int i=0; i<n; i++ )
			{
				y.push_back(i+M_PI);
				z.push_back(i+cos(1));
			}
			std::vector< std::vector<double> > A;
			A.push_back(x);
			A.push_back(y);
			A.push_back(z);
			for( int i=0; i<3; i++ )
			{
				for( int j=0; j<n; j++ )
				{
					std::cout << "A[" << i << "][" << j << "] = " << A[i][j] << std::endl;
				}
			}
			std::vector< std::vector<std::vector<double>> > B;
			B.push_back(A);
			for( int i=0; i<3; i++ )
			{
				for( int j=0; j<n; j++ )
				{
					std::cout << "B[0][" << i << "][" << j << "] = " << B[0][i][j] << std::endl;
				}
			}
		}
	}
	return 0;
}
