//mainScript.cpp
#include <iostream>
#include <string>
#include "fromFile.hpp"

int main()
{
	int d, m, n;
	std::string fileName( "in.txt" );
	d = getDims( fileName );
	m = getRows( fileName );
	if( d == 2 )
	{
		n = getCols( fileName );
		double x[m*n];
		getMat( fileName, x );
		std::cout << "d = " << d << std::endl;
		std::cout << "m = " << m << std::endl;
		std::cout << "n = " << n << std::endl;
		for( int i=0; i<m*n; i++ )
		{
			std::cout << "x[" << i << "] = " << x[i] << std::endl;
		}
	}
	else if( d == 1 )
	{
		double x[m];
		getVec( fileName, x );
		std::cout << "d = " << d << std::endl;
		std::cout << "m = " << m << std::endl;
		for( int i=0; i<m; i++ )
		{
			std::cout << "x[" << i << "] = " << x[i] << std::endl;
		}
	}
	else
	{
		std::cerr << "Error:  Only supporting d=1 and d=2." <<std::endl;
		return(EXIT_FAILURE);
	}
}
