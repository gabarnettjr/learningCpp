#ifndef CAR_HPP
#define CAR_HPP

#include <iostream>

namespace CPPbook{

class Car{

    protected:

        int km;

    public:

        Car( int d = 0 ) : km(d) {
        }

        virtual void travel( int d ) {
            km += d;
        }

        virtual void printTraveled() {
            std::cout << "The car has traveled " << km << " km." << std::endl;
        }

        virtual ~Car() {
        }

};
} //namespace CPPbook

#endif //CAR_HPP
