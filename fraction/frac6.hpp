//frac6.hpp
#ifndef FRACTION_HPP
#define FRACTION_HPP

#include <iostream>

namespace CPPbook
{
	class Fraction
	{
		private:

			int numer;

			int denom;

		public:

			Fraction( int=0, int=1 );

			friend Fraction operator * ( const Fraction&, const Fraction& );

			const Fraction& operator *= ( const Fraction& );

			friend bool operator < ( const Fraction&, const Fraction& );

			void printOn( std::ostream& );

			void scanFrom( std::istream& );
	};

	//inline Fraction operator * ( const Fraction& a, const Fraction& b )
	//{
	//	return Fraction( a.numer * b.numer, a.denom * b.denom );
	//}

	inline bool operator < (const Fraction& a, const Fraction& b )
	{
		return a.numer * b.denom < b.numer * a.denom;
	}

	inline std::ostream& operator << ( std::ostream& strm, const Fraction& f )
	{
		f.printOn( strm );
		return strm;
	}

	inline std::istream& operator >> ( std::istream& strm, Fraction& f )
	{
		f.scanFrom( strm );
		return strm;
	}
}

#endif
