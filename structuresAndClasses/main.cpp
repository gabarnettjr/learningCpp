#include <iostream>
#include "metal.hpp"
#include "point.hpp"

metal changeName( metal a ) {
    a.name = "hahaha";
    return a;
}

int main()
{
    metal a;
    a.name = "chromium";
    a.density = 7.19;
    a.meltPt = 1890;
    a.tensileMod = 289;
    a.daysDeliv = 19;

    a.printInfo();
    a = changeName( a );
    a.printInfo();

    const Point p( 1.6, -3.23 );
    const double x = p.getX();
    const double y = p.getY();
    p.prnt();
    
    Point q;
    double z = q.getX();
    double w = q.getY();
    q.prnt();

    return 0;
}
