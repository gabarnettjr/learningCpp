#include "cross.hpp"

int crossSum( long number )
{
	int cross = 0;

	while( number > 0 )
	{
		cross = cross + number % 10;
		number = number / 10;
	}
	return cross;
}
