//fromFile.cpp
#include <iostream>
#include <fstream>
#include <string>

void getVec( std::string fileName, double x[] )
{
	std::ifstream inFile;
	inFile.open( fileName );
	int d, m;
	inFile >> d >> m;
	for( int i=0; i<m; i++ )
	{
		inFile >> x[i];
	}
	inFile.close();
}

void getMat( std::string fileName, double x[] )
{
	std::ifstream inFile;
	inFile.open( fileName );
	int m, n;
	inFile >> m >> n;
	for( int i=0; i<m*n; i++ )
	{
		inFile >> x[i];
	}
	inFile.close();
}

void getArray( std::string fileName, int& d, int& m, int& n, double x[] )
{
	int N;
	std::ifstream inFile;
	inFile.open( fileName );
	inFile >> d >> m;
	if( d == 1 )
	{
		N = m;
	}
	else if( d == 2 )
	{
		inFile >> n;
		N = m*n;
	}
	for( int i=0; i<N; i++ )
	{
		inFile >> x[i];
	}
	inFile.close();
}

int getDims( std::string fileName )
{
	int d;
	std::ifstream inFile;
	inFile.open( fileName );
	inFile >> d;
	inFile.close();
	return d;
}

int getRows( std::string fileName )
{
	int d, m;
	std::ifstream inFile;
	inFile.open( fileName );
	inFile >> d >> m;
	inFile.close();
	return m;
}

int getCols( std::string fileName )
{
	int d, m, n;
	std::ifstream inFile;
	inFile.open( fileName );
	inFile >> d >> m >> n;
	inFile.close();
	return n;
}
