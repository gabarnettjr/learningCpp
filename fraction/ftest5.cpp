#include <iostream>
#include <cstdlib>
#include "frac5.hpp"

int main( int argc, char* argv[] )
{
	using namespace CPPbook;
	
	const Fraction a(7,3);
	Fraction x;
	
	std::cout << a << std::endl;
	
	std::cout << "Enter fraction (numer/denom):" << std::endl;
	if( ! ( std::cin >> x ) )
	{
		std::cerr << "Error during input of fraction." << std::endl;
		return EXIT_FAILURE;
	}
	std::cout << "Input was: " << x << std::endl;
	
	while( x <= 1000 )
	{
		x =  x*a;
		std::cout << x << std::endl;
	}
}
