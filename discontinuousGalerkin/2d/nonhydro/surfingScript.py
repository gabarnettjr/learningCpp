import numpy
import matplotlib.pyplot as plt

np = int( numpy.loadtxt( "np.txt" ) )
ne = int( numpy.loadtxt( "ne.txt" ) )
nTimesteps = int( numpy.loadtxt("nTimesteps.txt") )
t = float( numpy.loadtxt( "t.txt" ) )
dt = float( numpy.loadtxt( "dt.txt" ) )
x = numpy.loadtxt( "x.txt" )
z = numpy.loadtxt( "z.txt" )

a = min(x)
b = max(x)
dz = z[1] - z[0]
c = min(z)-dz/2.
d = max(z)+dz/2.
nLev = len(z)
n = np * ne
delta = 20                                           #number of timesteps between plots
zz = numpy.transpose( numpy.tile( z, (np,1) ) )

for i in numpy.arange(0,nTimesteps+1,delta) :
    rho = numpy.loadtxt( './snapshots/' + str(i).zfill(6) + '.txt' )
    if i != 0 :
        t = t + delta*dt
    if i == nTimesteps :
        for j in range(ne) :
            xx = numpy.tile( x[j*np:(j+1)*np], (nLev,1) ) 
            approx = rho[j*np:(j+1)*np]
            for k in numpy.arange(1,nLev) :
                approx = numpy.vstack(( approx, rho[k*n+j*np:k*n+(j+1)*np] ))
            exact = numpy.exp( -10 * ( xx**2 + zz**2 ) )
            plt.contourf( xx, zz, approx-exact, numpy.linspace(-.06,.06,13) )
        plt.axis( [a,b,c,d] )
        plt.title( 'approx - exact, t = {0:02.3f}'.format(t) )
        plt.colorbar()
        plt.waitforbuttonpress()
        plt.clf()
    else :
        for j in range(ne) :
            xx = numpy.tile( x[j*np:(j+1)*np], (nLev,1) ) 
            approx = rho[j*np:(j+1)*np]
            for k in numpy.arange(1,nLev) :
                approx = numpy.vstack(( approx, rho[k*n+j*np:k*n+(j+1)*np] ))
            plt.contourf( xx, zz, approx, numpy.arange(-.15,1.25,.1) )
        plt.axis( [a,b,c,d] )
        plt.title( 'numerical solution, t = {0:02.3f}'.format(t) )
        plt.colorbar()
        plt.waitforbuttonpress()
        plt.clf()