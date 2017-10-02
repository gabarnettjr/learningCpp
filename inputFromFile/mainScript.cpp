//mainScript.cpp
#include <iostream>
#include <string>
#include <fstream>
#include "fromFile.hpp"

int main()
{
	std::string fileName( "in.txt" );
	std::ifstream inFile;
	inFile.open( fileName );
	int m, n;
	inFile >> m >> n;
	inFile.close();
	double x[m*n];
	getMat( fileName, x );
	
	std::cout << "m = " << m << std::endl;
	std::cout << "n = " << n << std::endl;
	for( int i=0; i<m*n; i++ )
	{
		std::cout << "x[" << i << "] = " << x[i] << std::endl;
	}
}
