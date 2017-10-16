#include <iostream>
#include "cross.hpp"

void printCrossSum( long );

int main()
{
	printCrossSum( 12345678 );

	printCrossSum( 0 );

	printCrossSum( 13*77 );
}

void printCrossSum( long number )
{
	std::cout << "The cross sum of " << number << " is " << crossSum(number) << std::endl;
}
