import numpy
import matplotlib.pyplot as plt

nTimesteps = int( numpy.loadtxt("nTimesteps.txt") )
t = float( numpy.loadtxt( "t.txt" ) )
dt = float( numpy.loadtxt( "dt.txt" ) )

Rd = float( numpy.loadtxt( "Rd.txt" ) )
Po = float( numpy.loadtxt( "Po.txt" ) )
Cp = float( numpy.loadtxt( "Cp.txt" ) )
Cv = float( numpy.loadtxt( "Cv.txt" ) )

with open("s.txt") as testCase :
    s = testCase.readlines();
    s = s[0];

x = numpy.loadtxt( "x.txt" )
z = numpy.loadtxt( "z.txt" )
rhoBar = numpy.loadtxt( "rhoBar.txt" )
pBar = numpy.loadtxt( "pBar.txt" )
thetaBar = numpy.loadtxt( "thetaBar.txt" )
piBar = numpy.loadtxt( "piBar.txt" )
saveDelta = int( numpy.loadtxt( "saveDelta.txt" ) )    #number of timesteps between saves

###########################################################################

plotDelta = 5 * saveDelta                              #number of timesteps between plots
a = min(x)
b = max(x)
dx = x[1] - x[0]
dz = z[1] - z[0]
c = min(z) - dz/2.
d = max(z) + dz/2.
n = len(x)
nLev = len(z)
xx = numpy.tile( x, (nLev,1) )
zz = numpy.transpose( numpy.tile( z, (n,1) ) )

var = "thetaPrime"                                     #choose which variable to look at

fig = plt.figure()

for i in numpy.arange(0,nTimesteps+1,plotDelta) :

    if var == "thetaPrime" :
        rho = numpy.loadtxt( './rho/' + str(i).zfill(6) + '.txt' )
        tmp = numpy.loadtxt( './rhoTh/' + str(i).zfill(6) + '.txt' )
        tmp = tmp / rho - thetaBar
        if s == "risingBubble" :
            cv = numpy.arange( -.15, 2.25, .1 )
        elif s == "inertiaGravityWaves" :
            cv = numpy.arange( -.0015, .0035, .0005 )
        elif ( s == "densityCurrent" ) | ( s == "movingDensityCurrent" ) :
            cv = numpy.arange( -16.5, 1.5, 1 )
    elif var == "rhoPrime" :
        rho = numpy.loadtxt( './rho/' + str(i).zfill(6) + '.txt' )
        tmp = rho - rhoBar
        if s == "risingBubble" :
            cv = numpy.arange( -.021, .002, .001 );
    elif var == "pPrime" :
        rhoTh = numpy.loadtxt( './rhoTh/' + str(i).zfill(6) + '.txt' )
        tmp = Po * ( Rd * rhoTh / Po ) ** (Cp/Cv) - pBar
        if s == "risingBubble" :
            cv = numpy.arange( -82.5, 87.5, 5 )
        elif s == "densityCurrent" :
            cv = numpy.arange( -525, 575, 50 )
        elif s == "inertiaGravityWaves" :
            cv = numpy.arange( -10, 10, 1 )
    elif var == "piPrime" :
        tmp = numpy.loadtxt( './rhoTh/' + str(i).zfill(6) + '.txt' )
        tmp = Po * ( Rd * tmp / Po ) ** (Cp/Cv)
        tmp = ( tmp / Po ) ** (Rd/Cp) - piBar
        if s == "risingBubble" :
            cv = numpy.arange( -.0081, -.0059, .0001 )
    elif var == "u" :
        rho = numpy.loadtxt( './rho/' + str(i).zfill(6) + '.txt' )
        tmp = numpy.loadtxt( './rhoU/' + str(i).zfill(6) + '.txt' )
        tmp = tmp / rho
        if s == "risingBubble" :
            cv = numpy.arange( -6.5, 7.5, 1 )
        elif s == "inertiaGravityWaves" :
            cv = numpy.arange( 19.99, 20.011, .001 )
    elif var == "w" :
        rho = numpy.loadtxt( './rho/' + str(i).zfill(6) + '.txt' )
        tmp = numpy.loadtxt( './rhoW/' + str(i).zfill(6) + '.txt' )
        tmp = tmp / rho
        if s == "risingBubble" :
            cv = numpy.arange( -8.5, 11.5, 1 )
        elif s == "inertiaGravityWaves" :
            cv = numpy.arange( -.01, .0101, .001 )
    
    tmp = numpy.reshape( tmp, (nLev,n) )
    
    if i != 0 :
        t = t + plotDelta*dt
    plt.contourf( xx, zz, tmp, cv )
    
    plt.axis( [a,b,c,d] )
    plt.title( 'numerical solution, t = {0:02.0f}'.format(t) )
    plt.colorbar()
    fig.savefig( './pngs/'+'{0:04d}'.format(numpy.int(numpy.round(t)+1e-12))+'.png', bbox_inches = 'tight' )
    plt.clf()
