#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "Parameters.hpp"
#include "FDmesh.hpp"
#include "Variables.hpp"
#include "TimeStepper.hpp"

//Numerical solution for the 2D nonhydrostatic governing equations using
//finite differences laterally and finite differences vertically.
//Currently all variables are collocated at layer midpoints.

//Instructions:
// (*) mkdir u w dpids th
// (*) g++ Parameters.hpp FDmesh.hpp Variables.hpp TimeStepper.hpp main.cpp
// (*) ./a
// (*) python surfingScript.py

int main()
{
    //choose which test case to do with this string (see Parameters.hpp for choices)
    const std::string st( "risingBubble" );
    
    //Initialize parameters based on test case:
    const Parameters P( st );
    
    //Initialize mesh based on parameters:
    const FDmesh M( st, P );
    
    //Initialize variables based on test case, parameters, and mesh:
    Variables V( st, P, M );
    
    //set Runge-Kutta time stepper:
    TimeStepper T( P, V );
    
    int i, j, k;
    
    ///////////////////////////////////////////////////////////////////////
    
    //Open output file stream for saving things:
    std::ofstream outFile;
    
    //Precision for printing to text files:
    const int pr = 16;
    
    //save name of test:
    outFile.open( "s.txt" );
    outFile << s;
    outFile.close();
    
    //Save start time:
    outFile.open( "t.txt" );
    outFile << std::scientific << std::setprecision(pr) << P.t;
    outFile.close();
    
    //Save time increment:
    outFile.open( "dt.txt" );
    outFile << P.dt;
    outFile.close();
    
    //Save number of time-steps:
    outFile.open( "nTimesteps.txt" );
    outFile << P.nTimesteps;
    outFile.close();
    
    //save number of timesteps between saves:
    outFile.open( "saveDelta.txt" );
    outFile << P.saveDelta;
    outFile.close();
    
    //save array of x-coordinates:
    outFile.open( "x.txt" );
    for( i = 0; i < P.n; i++ ) {
        outFile << M.x[i] << " ";
    }
    outFile.close();
    
    //save array of z-coordinates:
    outFile.open( "z.txt" );
    for( i = 0; i < P.nLev; i++ ) {
        outFile << M.z[i] << " ";
    }
    outFile.close();
    
    //save background density:
    outFile.open( "rhoBar.txt" );
    for( i = 0; i < P.N; i++ ) {
        outFile << V.rhoBar[i] << " ";
    }
    outFile.close();
    
    //save background potential temperature"
    outFile.open( "thetaBar.txt" );
    for( i = 0; i < P.N; i++ ) {
        outFile << V.thetaBar[i] << " ";
    }
    outFile.close();
    
    //save background Exner pressure:
    outFile.open( "piBar.txt" );
    for( i = 0; i < P.N; i++ ) {
        outFile << V.piBar[i] << " ";
    }
    outFile.close();
    
    //save background pressure:
    outFile.open( "pBar.txt" );
    for( i = 0; i < P.N; i++ ) {
        outFile << V.pBar[i] << " ";
    }
    outFile.close();
    
    //save constants:
    outFile.open( "Rd.txt" );
    outFile << P.Rd;
    outFile.close();
    outFile.open( "Po.txt" );
    outFile << P.Po;
    outFile.close();
    outFile.open( "Cp.txt" );
    outFile << P.Cp;
    outFile.close();
    outFile.open( "Cv.txt" );
    outFile << P.Cv;
    outFile.close();
    
    ///////////////////////////////////////////////////////////////////////
    
    //Time stepping:
    
    // double rhoMin;
    // double rhoMax;
    double pMin;
    double pMax;
    // double pipMin;
    // double pipMax;
    double thpMin;
    double thpMax;
    double uMin;
    double uMax;
    double wMin;
    double wMax;
    
    for( j = 0; j < P.nTimesteps+1; j++ ) {
        
        if( j%P.saveDelta == 0 ) {
        
            //save rho:
            std::stringstream s_rho;
            s_rho << "./rho/" << std::setfill('0') << std::setw(6) << j << ".txt";
            outFile.open( s_rho.str() );
            for( i = 0; i < P.N; i++ ) {
                outFile << V.rho[i] << " ";
            }
            outFile.close();
            //save rhoU:
            std::stringstream s_rhoU;
            s_rhoU << "./rhoU/" << std::setfill('0') << std::setw(6) << j << ".txt";
            outFile.open( s_rhoU.str() );
            for( i = 0; i < P.N; i++ ) {
                outFile << V.rhoU[i] << " ";
            }
            outFile.close();
            //save rhoW:
            std::stringstream s_rhoW;
            s_rhoW << "./rhoW/" << std::setfill('0') << std::setw(6) << j << ".txt";
            outFile.open( s_rhoW.str() );
            for( i = 0; i < P.N; i++ ) {
                outFile << V.rhoW[i] << " ";
            }
            outFile.close();
            //save rhoTh:
            std::stringstream s_rhoTh;
            s_rhoTh << "./rhoTh/" << std::setfill('0') << std::setw(6) << j << ".txt";
            outFile.open( s_rhoTh.str() );
            for( i = 0; i < P.N; i++ ) {
                outFile << V.rhoTh[i] << " ";
            }
            outFile.close();
        
            //print min and max:
            std::cout << "t = " << V.t << std::endl;
            // pipMin = pow( V.p[0]/P.Po, P.Rd/P.Cp ) - V.piBar[0];
            // pipMax = pow( V.p[0]/P.Po, P.Rd/P.Cp ) - V.piBar[0];
            pMin = V.p[0] - V.pBar[0];
            pMax = V.p[0] - V.pBar[0];
            // rhoMin = V.rho[0] - V.rhoBar[0];
            // rhoMax = V.rho[0] - V.rhoBar[0];
            thpMin = V.rhoTh[0] / V.rho[0] - V.thetaBar[0];
            thpMax = V.rhoTh[0] / V.rho[0] - V.thetaBar[0];
            uMin = V.rhoU[0] / V.rho[0];
            uMax = V.rhoU[0] / V.rho[0];
            wMin = V.rhoW[0] / V.rho[0];
            wMax = V.rhoW[0] / V.rho[0];
            for( k = 1; k < P.N; k++ ) {
                // if( pow(V.p[k]/P.Po,P.Rd/P.Cp)-V.piBar[k] < pipMin ) {
                    // pipMin = pow(V.p[k]/P.Po,P.Rd/P.Cp)-V.piBar[k];
                // }
                // if( pow(V.p[k]/P.Po,P.Rd/P.Cp)-V.piBar[k] > pipMax ) {
                    // pipMax = pow(V.p[k]/P.Po,P.Rd/P.Cp)-V.piBar[k];
                // }
                if( V.p[k]-V.pBar[k] < pMin ) {
                    pMin = V.p[k]-V.pBar[k];
                }
                if( V.p[k]-V.pBar[k] > pMax ) {
                    pMax = V.p[k]-V.pBar[k];
                }
                // if( V.rho[k]-V.rhoBar[k] < rhoMin ) {
                    // rhoMin = V.rho[k]-V.rhoBar[k];
                // }
                // if( V.rho[k]-V.rhoBar[k] > rhoMax ) {
                    // rhoMax = V.rho[k]-V.rhoBar[k];
                // }
                if( V.rhoTh[k]/V.rho[k]-V.thetaBar[k] < thpMin ) {
                    thpMin = V.rhoTh[k]/V.rho[k]-V.thetaBar[k];
                }
                if( V.rhoTh[k]/V.rho[k]-V.thetaBar[k] > thpMax ) {
                    thpMax = V.rhoTh[k]/V.rho[k]-V.thetaBar[k];
                }
                if( V.rhoU[k]/V.rho[k] < uMin ) {
                    uMin = V.rhoU[k]/V.rho[k];
                }
                if( V.rhoU[k]/V.rho[k] > uMax ) {
                    uMax = V.rhoU[k]/V.rho[k];
                }
                if( V.rhoW[k]/V.rho[k] < wMin ) {
                    wMin = V.rhoW[k]/V.rho[k];
                }
                if( V.rhoW[k]/V.rho[k] > wMax ) {
                    wMax = V.rhoW[k]/V.rho[k];
                }
            }
            // std::cout << "minPip = " << pipMin << std::endl;
            // std::cout << "maxPip = " << pipMax << std::endl;
            std::cout << "minP = " << pMin << std::endl;
            std::cout << "maxP = " << pMax << std::endl;
            // std::cout << "minRho = " << rhoMin << std::endl;
            // std::cout << "maxRho = " << rhoMax << std::endl;
            std::cout << "minThp = " << thpMin << std::endl;
            std::cout << "maxThp = " << thpMax << std::endl;
            std::cout << "minU = " << uMin << std::endl;
            std::cout << "maxU = " << uMax << std::endl;
            std::cout << "minW = " << wMin << std::endl;
            std::cout << "maxW = " << wMax << std::endl << std::endl;
        }
        
        //advance (rho,rhoU,rhoW,rhoTh,P) and t with a single Runge-Kutta time step:
        T.rk( P, V );
    }
    
    return 0;
}
