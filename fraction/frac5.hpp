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
		
			Fraction();
			Fraction( int );
			Fraction( int, int );
			// Fraction( int = 0, int = 1 );	//weird default thing
			
			Fraction operator * ( const Fraction& ) const;
			
			const Fraction& operator *= ( const Fraction& );
			
			bool operator < ( const Fraction& ) const;
			
			bool operator <= ( const Fraction& ) const;
			
			bool operator > ( const Fraction& ) const;
			
			bool operator >= ( const Fraction& ) const;
			
			void printOn( std::ostream& ) const;
			
			void scanFrom( std::istream& );
			
	};
	
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