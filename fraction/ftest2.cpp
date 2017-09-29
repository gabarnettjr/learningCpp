#include "frac.hpp"

int main()
{
	CPPbook::Fraction x;
	CPPbook::Fraction w(7,3);
	
	x = w * w;
	
	while( x < CPPbook::Fraction(1000) )
	{
		x *= w;
		x.print();
	}
	
	CPPbook::Fraction z(x>w);
	
	w.print();
	z.print();
}