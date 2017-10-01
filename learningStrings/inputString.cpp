//Learning how to use information in the separate
//file, in.txt, to create an array.
//The first two numbers in in.txt give the size
//of the matrix, and the rest of the numbers
//give the entries in the matrix, so the first
//two need to be ints and the remaining can
//be doubles.

#include <iostream>
#include <fstream>

int main()
{
	int d;
	std::ifstream inFile;
	inFile.open( "in.txt" );
	inFile >> d;
	if( d == 1 )
	{
		int m;
		inFile >> m;
		std::cout << "d = " << d << std::endl;
		std::cout << "m = " << m << std::endl;
		double x[m];
		for( int i=0; i<m; i++ )
		{
			inFile >> x[i];
		}
		inFile.close();
		for( int i=0; i<m; i++ )
		{
			std::cout << "x[" << i << "] = " << x[i] << std::endl;
		}
	}
	else if( d == 2 )
	{
		int m, n;
		inFile >> m >> n;
		std::cout << "d = " << d << std::endl;
		std::cout << "m = " << m << std::endl;
		std::cout << "n = " << n << std::endl;
		double x[m*n];
		for( int i=0; i<m*n; i++ )
		{
			inFile >> x[i];
		}
		inFile.close();
		for( int i=0; i<m*n; i++ )
		{
			std::cout << "x[" << i << "] = " << x[i] << std::endl;
		}
	}
	//double x[m*n];
	//for( int i=0; i<m*n; i++ )
	//{
	//	inFile >> x[i];
	//}
	//inFile.close();
	//for( int i=0; i<m*n; i++ )
	//{
	//	std::cout << "x[" << i << "] = " << x[i] << std::endl;
	//}
}
