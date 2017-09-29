#include <iostream>
#include <string>
#include <cstdlib>

int main( int argc, char* argv[] )
{
	std::string progName = argv[0];
	
	if( argc > 1 )
	{
		std::cout << progName << std::endl << "has " << argc-1 << " parameters:" << std::endl;
		double x = strtod( argv[1], NULL );
		double y = strtod( argv[2], NULL );
		std::cout << "argv[1] + argv[2] = " << x+y << std::endl;
	}
	else
	{
		std::cout << progName << std::endl << "was called without parameters." << std::endl;
	}
	for( int i=0; i<argc; ++i )
	{
		std::cout << "argv[" << i << "]: " << argv[i] << std::endl;
	}
}
