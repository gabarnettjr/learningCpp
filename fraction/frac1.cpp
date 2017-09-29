#include <iostream>
#include <cstdlib>
#include "frac.hpp"

namespace CPPbook
{

	Fraction::Fraction()
		: numer(0), denom(1)
	{
	}

	Fraction::Fraction( int n )
		: numer(n), denom(1)
	{
	}

	Fraction::Fraction( int n, int d )
		: numer(n), denom(d)
	{
		if( d == 0 )
		{
			std::cerr << "error: denominator is 0" << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

	void Fraction::print()
	{
		std::cout << numer << '/' << denom << std::endl;
	}

}