#ifndef COMPLEX_HPP
#define COMPLEX_HPP

class Complex {

    public:

        Complex();
        Complex( double );
        Complex( double, double );

        double getRe() const;
        double getIm() const;

        Complex operator+( const Complex& ) const;
        Complex operator-( const Complex& ) const;
        Complex operator*( const Complex& ) const;
        Complex& operator+=( const Complex& );
        Complex& operator-=( const Complex& );

    private:

        double re;
        double im;

};

Complex::Complex() {
    re = 0.;
    im = 0.;
}

Complex::Complex( double a ) {
    re = a;
    im = 0.;
}

Complex::Complex( double a, double b ) {
    re = a;
    im = b;
}

double Complex::getRe() const {
    return re;
}

double Complex::getIm() const {
    return im;
}

Complex Complex::operator+( const Complex& c ) const {
    Complex ans( re + c.getRe(), im + c.getIm() );
    return ans;
}

Complex Complex::operator-( const Complex& c ) const {
    Complex ans( re - c.getRe(), im - c.getIm() );
    return ans;
}

Complex Complex::operator*( const Complex& c ) const {
    Complex ans( re*c.getRe() - im*c.getIm(), re*c.getIm() + im*c.getRe() );
    return ans;
}

Complex& Complex::operator+=( const Complex& c ) {
    re = re + c.getRe();
    im = im + c.getIm();
    return *this;       //optional?
}

Complex& Complex::operator-=( const Complex& c ) {
    re = re - c.getRe();
    im = im - c.getIm();
    return *this;       //optional?
}

#endif
