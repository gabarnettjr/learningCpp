#include <iostream>
#include <cstdlib>
#include "frac.hpp"

namespace CPPbook
{
	
	//Default constructor:
	Fraction::Fraction()
		: numer(0), denom(1)
	{
	}
	
	//Integer constructor:
	Fraction::Fraction( int n )
		: numer(n), denom(1)
	{
	}
	
	//Full fraction constructor:
	Fraction::Fraction( int n, int d )
		: numer(n), denom(d)
	{
		if( d == 0 )
		{
			std::cerr << "error: denominator is 0" << std::endl;
			std::exit(EXIT_FAILURE);
		}
		if ( d < 0 )
		{
			numer = -n;
			denom = -d;
		}
		else
		{
			numer = n;
			denom = d;
		}
	}
	
	//Print fraction:
	void Fraction::print()
	{
		std::cout << numer << '/' << denom << std::endl;
	}
	
	//multiply by fraction:
	Fraction Fraction::operator * (Fraction f )
	{
		return Fraction( numer * f.numer, denom * f.denom );
	}
	
	//set fraction equal to itself times incoming fraction:
	Fraction Fraction::operator *= (Fraction f )
	{
		*this = *this * f;
		return *this;
	}
	
	//check if fraction is less than incoming fraction:
	bool Fraction::operator < (Fraction f )
	{
		return numer*f.denom < f.numer*denom;
	}
	
	//check if fraction is greater than incoming fraction:
	bool Fraction::operator > (Fraction f )
	{
		return numer*f.denom > f.numer*denom;
	}

}