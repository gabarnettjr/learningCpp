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
		
			Fraction( int = 0, int = 1 );
			
			void print() const;
			
			Fraction operator * ( const Fraction& ) const;
			
			const Fraction& operator *= ( const Fraction& );
			
			bool operator < ( const Fraction& ) const;
			
			bool operator > ( const Fraction& ) const;
			
	};
	
	// inline void Fraction::print() const
	// {
		// std::cout << numer << '/' << denom << std::endl;
	// }
	
	// inline Fraction Fraction::operator * ( const Fraction& f ) const
	// {
		// return Fraction( numer*f.numer, denom*f.denom );
	// }
	
}

#endif