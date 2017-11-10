#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "Constants.hpp"
#include "DGmesh.hpp"
#include "Variables.hpp"
#include "TimeStepper.hpp"

//Numerical solution for the 2D nonhydrostatic governing equations using
//discontinuous Galerkin laterally and finite differences vertically.
//Currently all variables are collocated at layer midpoints.

//These are the equations being solved:
//  (1) d(rho)/dt    = -d(rho*u)/dx - d(rho*w)/dz,
//  (2) d(rho*u)/dt  = -d(rho*u*u+P)/dx - d(rho*u*w)/dz,
//  (3) d(rho*w)/dt  = -d(rho*w*u)/dx - d(rho*w*w+P)/dz - rho*g,
//  (4) d(rho*th)/dt = -d(rho*th*u)/dx - d(rho*th*w)/dz,
//  (5) P = Po*(rho*Rd*th/Po)^(Cp/Cv),
//with Po=10^5, Rd=287.04, Cp=1005, Cv=Cp-Rd, g=9.81.
//w=0 ( => dP/dz=-rho*g ) is enforced on the top and bottom boundaries.

//Instructions:
// (*) mkdir rho rhoU rhoW rhoTh
// (*) g++ Constants.hpp DGmesh.hpp Variables.hpp TimeStepper.hpp main.cpp
// (*) ./a
// (*) python surfingScript.py (NOT WORKING YET)

int main()
{   
    //Parameters that the user chooses:
    const double a = 0.;                          //left endpoint
    const double b = 10000.;                      //right endpoint
    const double c = 0.;                          //bottom endpoint
    const double d = 10000.;                      //top endpoint
    const int np = 4;                             //number of polynomials per element (2, 3, or 4)
    const int ne = 25;                            //number of elements per level (layer)
    const int nLev = 100;                         //number of levels (layers)
    const int rkStages = 4;                       //number of Runge-Kutta stages (2, 3, or 4)
    double t = 0.;                                //start time
    const double dt = 1./4.;                      //time increment
    const int nTimesteps = 6000;                 //number of timesteps
    const int saveDelta = 40;                     //number of time steps between saves
    
    int i, j, k;
    
    //physical constants:
    Constants C;
    
    //Initialize mesh:
    DGmesh M( a, b, c, d, np, ne, nLev );
    
    //set up the prognostic and diagnostic variables:
    Variables V( C, M );
    
    //set time stepper:
    TimeStepper T( C, M, V, rkStages, nTimesteps, t, dt );
    
    ///////////////////////////////////////////////////////////////////////
    
    //Open output file stream for saving things:
    std::ofstream outFile;
    
    //Precision for printing to text files:
    const int pr = 16;
    
    //Save number of polynomials per element:
    outFile.open( "np.txt" );
    outFile << np;
    outFile.close();
    
    //Save number of elements:
    outFile.open( "ne.txt" );
    outFile << ne;
    outFile.close();
    
    //Save start time:
    outFile.open( "t.txt" );
    outFile << std::scientific << std::setprecision(pr) << t;
    outFile.close();
    
    //Save time increment:
    outFile.open( "dt.txt" );
    outFile << dt;
    outFile.close();
    
    //Save number of time-steps:
    outFile.open( "nTimesteps.txt" );
    outFile << nTimesteps;
    outFile.close();
    
    //save number of timesteps between saves:
    outFile.open( "saveDelta.txt" );
    outFile << saveDelta;
    outFile.close();
    
    //save array of x-coordinates:
    outFile.open( "x.txt" );
    for( i = 0; i < M.n; i++ ) {
        outFile << M.x[i] << " ";
    }
    outFile.close();
    
    //save array of z-coordinates:
    outFile.open( "z.txt" );
    for( i = 0; i < nLev; i++ ) {
        outFile << M.z[i] << " ";
    }
    outFile.close();
    
    ///////////////////////////////////////////////////////////////////////
    
    //Time stepping:
    
    double pipMin;
    double pipMax;
    double thpMin;
    double thpMax;
    double uMin;
    double uMax;
    double wMin;
    double wMax;
    
    for( j = 0; j < nTimesteps+1; j++ ) {
        
        if( j%saveDelta == 0 ) {
        
            //save rho:
            std::stringstream s_rho;
            s_rho << "./rho/" << std::setfill('0') << std::setw(6) << j << ".txt";
            outFile.open( s_rho.str() );
            for( i = 0; i < M.N; i++ ) {
                outFile << V.rho[i] << " ";
            }
            outFile.close();
            //save rhoU:
            std::stringstream s_rhoU;
            s_rhoU << "./rhoU/" << std::setfill('0') << std::setw(6) << j << ".txt";
            outFile.open( s_rhoU.str() );
            for( i = 0; i < M.N; i++ ) {
                outFile << V.rhoU[i] << " ";
            }
            outFile.close();
            //save rhoW:
            std::stringstream s_rhoW;
            s_rhoW << "./rhoW/" << std::setfill('0') << std::setw(6) << j << ".txt";
            outFile.open( s_rhoW.str() );
            for( i = 0; i < M.N; i++ ) {
                outFile << V.rhoW[i] << " ";
            }
            outFile.close();
            //save rhoTh:
            std::stringstream s_rhoTh;
            s_rhoTh << "./rhoTh/" << std::setfill('0') << std::setw(6) << j << ".txt";
            outFile.open( s_rhoTh.str() );
            for( i = 0; i < M.N; i++ ) {
                outFile << V.rhoTh[i] << " ";
            }
            outFile.close();
        
            //print information:
            std::cout << "t = " << T.t << std::endl;
            pipMin = pow( V.P[0]/C.Po, C.Rd/C.Cp ) - V.piBar[0];
            pipMax = pow( V.P[0]/C.Po, C.Rd/C.Cp ) - V.piBar[0];
            thpMin = V.rhoTh[0] / V.rho[0] - V.thetaBar[0];
            thpMax = V.rhoTh[0] / V.rho[0] - V.thetaBar[0];
            uMin = V.rhoU[0] / V.rho[0];
            uMax = V.rhoU[0] / V.rho[0];
            wMin = V.rhoW[0] / V.rho[0];
            wMax = V.rhoW[0] / V.rho[0];
            for( k = 1; k < M.N; k++ ) {
                if( pow(V.P[k]/C.Po,C.Rd/C.Cp)-V.piBar[k] < pipMin ) {
                    pipMin = pow(V.P[k]/C.Po,C.Rd/C.Cp)-V.piBar[k];
                }
                if( pow(V.P[k]/C.Po,C.Rd/C.Cp)-V.piBar[k] > pipMax ) {
                    pipMax = pow(V.P[k]/C.Po,C.Rd/C.Cp)-V.piBar[k];
                }
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
            std::cout << "minPip = " << pipMin << std::endl;
            std::cout << "maxPip = " << pipMax << std::endl;
            std::cout << "minThp = " << thpMin << std::endl;
            std::cout << "maxThp = " << thpMax << std::endl;
            std::cout << "minU = " << uMin << std::endl;
            std::cout << "maxU = " << uMax << std::endl;
            std::cout << "minW = " << wMin << std::endl;
            std::cout << "maxW = " << wMax << std::endl << std::endl;
        }
        
        //advance (rho,rhoU,rhoW,rhoTh,P) and t with a single Runge-Kutta time step:
        T.rk( C, M, V );
    }
    
    return 0;
}
