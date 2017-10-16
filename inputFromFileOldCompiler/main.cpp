#include <iostream>
#include <string>
#include <fstream>

int main()
{
    //std::string mystring;
    //mystring = "This is a string";
    //std::cout << std::endl << mystring << std::endl;
    //
    //std::string name_inFile( "in.txt" );
    //std::string name_outFile( "out.txt" );
    char name_inFile[7] = "in.txt";
    char name_outFile[8] = "out.txt";

    std::ifstream inFile;
    inFile.open( name_inFile );
    int n;
    inFile >> n;
    //std::cout << std::endl << "n = " << n << std::endl;
    double x[n];
    for( int i=0; i<n; i++ ) {
        inFile >> x[i];
        //std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    inFile.close();

    std::ofstream outFile;
    outFile.open( name_outFile );
    outFile << n << std::endl << x[0]*x[0];
    for( int i=1; i<n; i++ ) {
        x[i] = x[i] * x[i];
        outFile << " " << x[i];
        //std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }
    outFile.close();

}
