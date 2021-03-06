import numpy
import matplotlib.pyplot as plt

np = int( numpy.loadtxt( "np.txt" ) )
ne = int( numpy.loadtxt( "ne.txt" ) )
nTimesteps = int( numpy.loadtxt("nTimesteps.txt") )
t = float( numpy.loadtxt( "t.txt" ) )
dt = float( numpy.loadtxt( "dt.txt" ) )

Rd = float( numpy.loadtxt( "Rd.txt" ) )
Po = float( numpy.loadtxt( "Po.txt" ) )
Cp = float( numpy.loadtxt( "Cp.txt" ) )
Cv = float( numpy.loadtxt( "Cv.txt" ) )

x = numpy.loadtxt( "x.txt" )
z = numpy.loadtxt( "z.txt" )
rhoBar = numpy.loadtxt( "rhoBar.txt" )
Pbar = numpy.loadtxt( "Pbar.txt" )
thetaBar = numpy.loadtxt( "thetaBar.txt" )
piBar = numpy.loadtxt( "piBar.txt" )
saveDelta = int( numpy.loadtxt( "saveDelta.txt" ) )    #number of timesteps between saves

###########################################################################

plotDelta = 10 * saveDelta;                            #number of timesteps between plots
a = min(x)
b = max(x)
dz = z[1] - z[0]
c = min(z)-dz/2.
d = max(z)+dz/2.
nLev = len(z)
n = np * ne
zz = numpy.transpose( numpy.tile( z, (np,1) ) )

for i in numpy.arange(0,nTimesteps+1,plotDelta) :

    # #thetaPrime:
    # rho = numpy.loadtxt( './rho/' + str(i).zfill(6) + '.txt' )
    # tmp = numpy.loadtxt( './rhoTh/' + str(i).zfill(6) + '.txt' )
    # tmp = tmp / rho - thetaBar
    # cv = numpy.arange( -.15, 2.25, .1 )
    
    # #rhoPrime:
    # rho = numpy.loadtxt( './rho/' + str(i).zfill(6) + '.txt' )
    # tmp = rho - rhoBar
    # cv = numpy.arange( -.021, .002, .001 );
    
    #Pprime:
    rhoTh = numpy.loadtxt( './rhoTh/' + str(i).zfill(6) + '.txt' )
    tmp = Po * ( Rd * rhoTh / Po ) ** (Cp/Cv) - Pbar
    cv = numpy.arange( -82.5, 87.5, 5 )
    
    # #piPrime:
    # tmp = numpy.loadtxt( './rhoTh/' + str(i).zfill(6) + '.txt' )
    # tmp = Po * ( Rd * tmp / Po ) ** (Cp/Cv)
    # tmp = ( tmp / Po ) ** (Rd/Cp) - piBar
    # cv = numpy.arange( -.0081, -.0059, .0001 )
    
    # #u:
    # rho = numpy.loadtxt( './rho/' + str(i).zfill(6) + '.txt' )
    # tmp = numpy.loadtxt( './rhoU/' + str(i).zfill(6) + '.txt' )
    # tmp = tmp / rho
    # cv = numpy.arange( -6.5, 7.5, 1 )
    
    # #w:
    # rho = numpy.loadtxt( './rho/' + str(i).zfill(6) + '.txt' )
    # tmp = numpy.loadtxt( './rhoW/' + str(i).zfill(6) + '.txt' )
    # tmp = tmp / rho
    # cv = numpy.arange( -8.5, 11.5, 1 )
    
    if i != 0 :
        t = t + plotDelta*dt
    for j in range(ne) :
        xx = numpy.tile( x[j*np:(j+1)*np], (nLev,1) ) 
        approx = tmp[j*np:(j+1)*np]
        for k in numpy.arange(1,nLev) :
            approx = numpy.vstack(( approx, tmp[k*n+j*np:k*n+(j+1)*np] ))
        plt.contourf( xx, zz, approx, cv )
    plt.axis( [a,b,c,d] )
    plt.title( 'numerical solution, t = {0:02.3f}'.format(t) )
    plt.colorbar()
    plt.waitforbuttonpress()
    plt.clf()