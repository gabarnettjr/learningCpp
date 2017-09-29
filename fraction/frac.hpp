#ifndef FRACTION_HPP
#define FRACTION_HPP

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
			void print();
			Fraction operator * (Fraction);
			Fraction operator *= (Fraction);
			bool operator < (Fraction);
			bool operator > (Fraction);
	};
}

#endif