#include "frac4.hpp"

int main()
{
	using namespace CPPbook;
	
	const Fraction w(7,3);
	
	w.print();
	
	Fraction x = w * w;
	
	while( x < Fraction(1000) )
	{
		x *= w;
		x.print();
	}
}