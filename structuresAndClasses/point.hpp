#ifndef POINT_HPP
#define POINT_HPP

#include <iostream>

class Point {

    public:

        Point();
        Point( const double, const double );

        const double getX() const;
        const double getY() const;

        void prnt() const;

    private:

        double x;
        double y;

};

Point::Point() {
    x = 0.;
    y = 0.;
}

Point::Point( const double a, const double b ) {
    x = a;
    y = b;
}

const double Point::getX() const {
    return x;
}

const double Point::getY() const {
    return y;
}

void Point::prnt() const {
    std::cout << std::endl;
    std::cout << "point = ( " << x << ", " << y << " )" << std::endl;
}

#endif
