#include <iostream>

int main()
{
	//initialize all integer variables:
	int x = 42;
	int y = 13;
	int z = 0;
	int w = 0;

	//long way to find the minimum of the two numbers:
	if( x < y )
	{
		z = x;
	}
	else
	{
		z = y;
	}

	//short way to find the minimum of the two numbers:
	w = x < y ? x : y;

	std::cout << "\nz = " << z << std::endl
		  << "w = "   << w << std::endl;
}
