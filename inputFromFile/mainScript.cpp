//mainScript.cpp
#include <iostream>
#include <string>
#include <fstream>
#include "fromFile.hpp"

int main()
{
	//read x in from "in.txt":
	std::string inputFileName( "in.txt" );
	int m, n;
	m = getRows( inputFileName );
	n = getCols( inputFileName );
	double x[m*n];
	getMat( inputFileName, x );
	
	//check and make sure data is correct:
	std::cout << "m = " << m << std::endl;
	std::cout << "n = " << n << std::endl;
	for( int i=0; i<m*n; i++ )
	{
		std::cout << "x[" << i << "] = " << x[i] << std::endl;
	}

	//replace x with -x:
	for( int i=0; i<m*n; i++ )
	{
		x[i] = -x[i];
	}

	//output negative array to "out.txt":
	std::ofstream outFile;
	outFile.open( "out.txt" );
	outFile << m << " " << n;
	for( int i=0; i<m*n; i++ )
	{
		outFile << " " << x[i];
	}
	outFile.close();
}
