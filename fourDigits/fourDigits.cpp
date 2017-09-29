#include <iostream>

int main()
{
	int counter = 0;
	int front = 0;
	int back = 0;
	int number = 0;

	for( number=1000; number<10000; ++number )
	{
		front = number / 100;
		back = number % 100;
		if( front*front + back*back == number )
		{
			std::cout << number << " = "
				  << front  << "*" << front << " + "
				  << back   << "*" << back << std::endl;
			++counter;
		}
	}
}
