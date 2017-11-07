import numpy as np
import matplotlib.pyplot as plt

t = float( np.loadtxt( "t.txt" ) )
dt = float( np.loadtxt( "dt.txt" ) )
x = np.loadtxt( "x.txt" )
z = np.loadtxt( "z.txt" )
dz = float( np.loadtxt( "dz.txt" ) )
a = min(x)
b = max(x)
c = min(z)-dz/2.
d = max(z)+dz/2.
nPoly = int( np.loadtxt( "np.txt" ) )
ne = int( np.loadtxt( "ne.txt" ) )
nLev = int( np.loadtxt( "nLev.txt" ) )
n = int( np.loadtxt( "n.txt" ) )
nTimesteps = int( np.loadtxt("nTimesteps.txt") )
delt = 20

zz = np.transpose( np.tile( z, (nPoly,1) ) )
#plt.ion()
for i in np.arange(0,nTimesteps+1,delt) :
    rho = np.loadtxt( './snapshots/' + str(i).zfill(6) + '.txt' )
    for j in range(ne) :
        xx = np.tile( x[j*nPoly:(j+1)*nPoly], (nLev,1) )
        tmp = rho[j*nPoly:(j+1)*nPoly];
        for k in np.arange(1,nLev) :
            tmp = np.vstack(( tmp, rho[k*n+j*nPoly:k*n+(j+1)*nPoly] ))
        plt.contourf( xx, zz, tmp )
    plt.axis( [a,b,c,d] )
    plt.title( '{0:02.3f}'.format(t) )
    t = t + delt*dt
    plt.waitforbuttonpress()
    plt.cla()

#plt.ioff()
#plt.plot( x, rho - np.exp(-10*x**2) )
#plt.plot( x, rho - np.cos(np.pi*x) )
#plt.plot( x, rho )
#plt.show()
