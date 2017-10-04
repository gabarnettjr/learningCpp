//frac6.hpp
#ifndef FRACTION_HPP
#define FRACTION_HPP

#include <iostream>

namespace CPPbook
{
	class Fraction
	{
		protected:

			int numer;

			int denom;

		public:

			class DenomIsZero
			{
			};
		
		//	Fraction();
		//	Fraction( int );
		//	Fraction( int, int );
			Fraction( int = 0, int = 1 );	//weird default thing
			
			friend Fraction operator * ( const Fraction& const Fraction& );
			
			const Fraction& operator *= ( const Fraction& );
			
			friend bool operator < ( const Fraction&, const Fraction& );
			
			friend bool operator <= ( const Fraction&, const Fraction& );
			
			friend bool operator > ( const Fraction&, const Fraction& );
			
			friend bool operator >= ( const Fraction&, const Fraction& );
			
			void printOn( std::ostream& ) const;
			
			void scanFrom( std::istream& );

			double toDouble() const;
	};

	inline Fraction operator * ( const Fraction& a, const Fraction& b )
	{
		return Fraction( a.numer * b.numer, a.denom * b.denom );
	}
	
	inline
	std::ostream& operator << ( std::ostream& strm, const Fraction& f )
	{
		f.printOn(strm);
		return strm;
	}
	
	inline
	std::istream& operator >> ( std::istream& strm, Fraction& f )
	{
		f.scanFrom(strm);
		return strm;
	}
	
}

#endif
