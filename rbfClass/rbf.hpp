#ifndef RBF_HPP
#define RBF_HPP

#include <iostream>
#include <cstdlib>

class rbf
{
	private:
		std::string rbfType;
		double parameter;
		double scaler;

	public:
	
		rbf();
		rbf( std::string, double );
		rbf( std::string, double, double );
		
		rbf operator d ( rbf& );
		
};

#endif