#include <iostream>

int main()
{
	int x = 0;
	do
	{
		++x;
		std::cout << "x: " << x << std::endl;
	}
	while( x < 7 );

	while( x < 15 )
	{
		std::cout << "x: " << x << std::endl;
		++x;
	}
}
