#include "frac.hpp"

int main()
{
	CPPbook::Fraction x;
	CPPbook::Fraction w(7,3);
	
	w.print();
	
	x = w;
	
	w = CPPbook::Fraction(100);
	
	x.print();
	w.print();
}