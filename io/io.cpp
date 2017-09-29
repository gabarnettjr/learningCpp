#include <iostream>
#include <cstdlib>

int main()
{
	int x, y;
	
	std::cout << "Integral division (x/y)\n\n";
	
	std::cout << "x: ";
	if( !( std::cin >> x ) )
	{
		std::cerr << "Error when reading an integer." << std::endl;
		return EXIT_FAILURE;
	}
	
	std::cout << "y: ";
	if( !( std::cin >> y ) )
	{
		std::cerr << "Error when reading an integer." << std::endl;
		return EXIT_FAILURE;
	}
	
	if( y == 0 )
	{
		std::cerr << "Error: division by zero." << std::endl;
		return EXIT_FAILURE;
	}
	
	std::cout << x << " divided by " << y << " gives " << x/y << "." << std::endl;
}