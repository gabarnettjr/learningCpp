#include <iostream>

int main()
{
	//for every character c with a value of 32 to 126:
	for ( unsigned char c=32; c<127; ++c )
	{
		//output value as a number and as a character:
		std::cout << "Value: "     << static_cast<int>(c)
			<< " Character: " << c
			<< std::endl;
	}
}
