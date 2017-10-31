#ifndef METAL_HPP
#define METAL_HPP

#include <iostream>
#include <string>

struct metal {
    std::string name;
    double density;
    int meltPt;
    int tensileMod;
    int daysDeliv;
    void printInfo() {
        std::cout << std::endl;
        std::cout << "name: " << name << std::endl;
        std::cout << "density: " << density << std::endl;
        std::cout << "melting point: " << meltPt << std::endl;
        std::cout << "tensile modulus: " << tensileMod << std::endl;
        std::cout << "days to delivery: " << daysDeliv << std::endl;
    }
};

#endif
