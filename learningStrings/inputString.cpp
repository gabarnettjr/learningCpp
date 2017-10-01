#include <iostream>
#include <fstream>

int main()
{
	int m;
	int n;
	std::ifstream inFile;
	inFile.open( "in.txt" );
	inFile >> m;
	inFile >> n;
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
