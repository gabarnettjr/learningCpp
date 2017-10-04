//rftest1.cpp
#include <iostream>
#include "rfrac1.hpp"

int main()
{
	CPPbook::RFraction x( 91, 39 );

	std::cout << x;
	std::cout << ( x.isReducible() ? " (reducible)" : " (non reducible)" ) std::endl;

	x.reduce();

	std::cout << x;
	std::cout << ( x.isReducible() ? " (reducible)" : " (non reducible)" ) std::endl;

	x *= 3;

	std::cout << x;
	std::cout << ( x.isReducible() ? " (reducible)" : " (non reducible)" ) std::endl;
}
