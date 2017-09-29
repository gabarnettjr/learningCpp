#include "frac4.hpp"

int main()
{
	using namespace CPPbook;
	
	const Fraction z(3,2);
	const Fraction w(7,3);
	
	std::cout << "z = ";
	z.print();
	
	std::cout << "w = ";
	w.print();
	
	const Fraction u = ( z >= w );
	u.print();
	
	Fraction x = w * w;
	
	while( x <= Fraction(1000) )
	{
		x *= w;
		x.print();
	}
}