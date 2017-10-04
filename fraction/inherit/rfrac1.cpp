//rfrac1.cpp
#include <algorithm>
#include <cstdlib>
#include "rfrac1.hpp"

namespace CPPbook
{
	unsigned RFraction::gcd() const
	{
		if( numer == 0 )
		{
			return denom;
		}
		unsigned divisor = std::min( std::abs(numer), denom );
		while( numer%divisor!=0 || denom%divisor!=0 )
		{
			--divisor;
		}
		return divisor;
	}

	void RFraction::reduce()
	{
		if( reducible )
		{
			int divisor = gcd();
			numer /= divisor;
			denom /= divisor;
			reducible = false;
		}
	}

	const RFraction& RFraction::operator *= ( const RFraction& f )
	{
		numer *= f.numer;
		denom *= f.denom;
		if( !reducible )
		{
			reducible = ( gcd() > 1 );
		}
		return *this;
	}

	void RFraction::scanFrom( std::istream& strm )
	{
		Fraction::scanFrom( strm );
		reducible = ( gcd() > 1 );
	}
}
