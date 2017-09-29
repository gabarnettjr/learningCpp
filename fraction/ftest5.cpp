#include "frac5.hpp"

int main()
{
	using namespace CPPbook;
	
	Fraction x;
	const Fraction z(3);
	const Fraction w(7,3);
	
	std::cout << "z = " << z << std::endl;
	
	std::cout << "w = " << w << std::endl;
	
	const Fraction u = ( z >= w );
	std::cout << "u = " << u << std::endl;
	
	x = w * w;
	
	while( x <= Fraction(1000) )
	{
		x =  x*w;
		std::cout << x << std::endl;
	}
}