import numpy
import matplotlib.pyplot as plt

np = int( numpy.loadtxt( "np.txt" ) )
ne = int( numpy.loadtxt( "ne.txt" ) )
nTimesteps = int( numpy.loadtxt("nTimesteps.txt") )
t = float( numpy.loadtxt( "t.txt" ) )
dt = float( numpy.loadtxt( "dt.txt" ) )
x = numpy.loadtxt( "x.txt" )
z = numpy.loadtxt( "z.txt" )
saveDelta = int( numpy.loadtxt( "saveDelta.txt" ) )    #number of timesteps between saves

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
    rho = numpy.loadtxt( './rho/' + str(i).zfill(6) + '.txt' )
    tmp = numpy.loadtxt( './rhoTh/' + str(i).zfill(6) + '.txt' )
    tmp = tmp / rho - 300.;
    if i != 0 :
        t = t + plotDelta*dt
    for j in range(ne) :
        xx = numpy.tile( x[j*np:(j+1)*np], (nLev,1) ) 
        approx = tmp[j*np:(j+1)*np]
        for k in numpy.arange(1,nLev) :
            approx = numpy.vstack(( approx, tmp[k*n+j*np:k*n+(j+1)*np] ))
        plt.contourf( xx, zz, approx, numpy.arange(-.15,2.15,.1) )
    plt.axis( [a,b,c,d] )
    plt.title( 'numerical solution, t = {0:02.3f}'.format(t) )
    plt.colorbar()
    plt.waitforbuttonpress()
    plt.clf()