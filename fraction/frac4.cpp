#include <cstdlib>
#include "frac4.hpp"

namespace CPPbook
{
	
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
	
	void Fraction::print() const
	{
		std::cout << numer << '/' << denom << std::endl;
	}
	
	Fraction Fraction::operator * ( const Fraction& f ) const
	{
		return Fraction( numer*f.numer, denom*f.denom );
	}
	
	//set fraction equal to itself times incoming fraction:
	const Fraction& Fraction::operator *= ( const Fraction& f )
	{
		*this = *this * f;
		return *this;
	}
	
	//check if fraction is less than incoming fraction:
	bool Fraction::operator < ( const Fraction& f ) const
	{
		return numer*f.denom < f.numer*denom;
	}
	
	bool Fraction::operator <= ( const Fraction& f ) const
	{
		return numer*f.denom <= f.numer*denom;
	}
	
	//check if fraction is greater than incoming fraction:
	bool Fraction::operator > ( const Fraction& f ) const
	{
		return numer*f.denom > f.numer*denom;
	}
	
	bool Fraction::operator >= ( const Fraction& f ) const
	{
		return numer*f.denom >= f.numer*denom;
	}

}