//frac6.cpp
#include <iostream>
#include <cstdlib>
#include "frac6.hpp"

namespace CPPbook
{
	
	////fraction constructor with no inputs:
	//Fraction::Fraction()
	//	: numer(0), denom(1)
	//{
	//}
	//
	////fraction constructor with one input:
	//Fraction::Fraction( int n )
	//	: numer(n), denom(1)
	//{
	//}
	
	//Full fraction constructor with two inputs:
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
	
	//function to multiply two fraction:
	Fraction Fraction::operator * ( const Fraction& a, const Fraction& b)
	{
		return Fraction( a.numer*b.numer, a.denom*b.denom );
	}
	
	//set fraction equal to itself times incoming fraction:
	const Fraction& Fraction::operator *= ( const Fraction& f )
	{
		*this = *this * f;
		return *this;
	}
	
	////check if fraction is less than incoming fraction:
	//bool Fraction::operator < ( const Fraction& f ) const
	//{
	//	return numer*f.denom < f.numer*denom;
	//}
	//
	////check if fraction is less than or equal to incoming fraction:
	//bool Fraction::operator <= ( const Fraction& f ) const
	//{
	//	return numer*f.denom <= f.numer*denom;
	//}
	//
	////check if fraction is greater than incoming fraction:
	//bool Fraction::operator > ( const Fraction& f ) const
	//{
	//	return numer*f.denom > f.numer*denom;
	//}
	//
	////check if fraction is greater than or equal to incoming fraction:
	//bool Fraction::operator >= ( const Fraction& f ) const
	//{
	//	return numer*f.denom >= f.numer*denom;
	//}
	
	void Fraction::printOn( std::ostream& strm ) const
	{
		strm << numer << '/' << denom;
	}
	
	void Fraction::scanFrom( std::istream& strm )
	{
		int n, d;
		strm >> n;
		if( strm.peek() == '/' )
		{
			strm.get();
			strm >> d;
		}
		else
		{
			d = 1;
		}
		if( ! strm )
		{
			return;
		}
		if( d == 0 )
		{
			strm.clear( strm.rdstate() | std::ios::failbit );
			return;
		}
		if( d < 0 )
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
	
}
