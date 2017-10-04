//ftest6.cpp
#include <iostream>
#include <cstdlib>
#include "frac6.hpp"

int main()
{
	const CPPbook::Fraction:: a( 7, 3 );
	CPPbook::Fraction x;


	std::cout << a << std::endl;

	std::cout << "enter fraction (numer/denom): ";
	if( !(std::cin>>x) )
	{
		std::cerr << "Error during input of fraction." << std::endl;
		return EXIT_FAILURE;
	}
	std::cout << "Input was: " << x << std::endl;

	while( x < 1000 )
	{
		x = x * a;
		std::cout << x << std::endl;
	}
}
