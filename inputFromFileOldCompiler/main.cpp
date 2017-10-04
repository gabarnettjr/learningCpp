#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

int main()
{
  std::string mystring;
  mystring = "This is a string";
  std::cout << mystring << std::endl;

  int m;
  //std::string name_inFile( "in.txt" );
  //std::string name_outFile( "out.txt" );
  char name_inFile[7] = "in.txt";
  char name_outFile[8] = "out.txt";

  std::ifstream inFile;
  inFile.open( name_inFile );
  inFile >> m;
  inFile.close();

  std::cout << "m = " << m << std::endl;

  int n = m*m;
  std::ofstream outFile;
  outFile.open( name_outFile );
  outFile << n;
  outFile.close();

  std::cout << "n = " << n << std::endl;
}
