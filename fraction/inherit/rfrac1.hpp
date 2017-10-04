//rfrac1.hpp
#ifndef RFRACTION_HPP
#define RFRACTION_HPP

#include "frac6.hpp"

namespace CPPbook
{
	class RFraction : public Fraction
	{
		protected:

			bool reducible;

			unsigned gcd() const;

		public:

			RFraction( int n = 0;  int d = 1 ) : Fraction( n, d )
			{
				reducible = ( gcd() > 1 );
			}

			const RFraction& operator *= ( const RFraction& );

			void scanFrom( std::istream& );

			void reduce();

			bool isReducible() const
			{
				return reducible;
			}
	}
}

#endif
