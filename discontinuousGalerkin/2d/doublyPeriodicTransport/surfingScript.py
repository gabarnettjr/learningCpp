import numpy as np
import matplotlib.pyplot as plt
import time

t = float( np.loadtxt("t.txt") )
dt = float( np.loadtxt("dt.txt") )
x = np.loadtxt("x.txt")
z = np.loadtxt("z.txt")
dz = np.loadtxt("dz.txt")
a = min(x)
b = max(x)
c = min(z)-dz/2;
d = max(z)+dz/2;
nTimesteps = int( np.loadtxt("nTimesteps.txt") )
delt = 20

plt.ion()
for i in np.arange(0,nTimesteps+1,delt) :
    rho = np.loadtxt('./snapshots/'+str(i).zfill(6)+'.txt')
    rho = np.reshape( rho, ( len(z), len(x) ) )
    plt.surf( xx, zz, rho )
    plt.axis( [a,b,-1.2,1.2] )
    plt.title( '{0:02.3f}'.format(t) )
    t = t + delt*dt
    plt.waitforbuttonpress()
    plt.cla()

plt.ioff()
plt.plot( x, rho - np.exp(-10*x**2) )
#plt.plot( x, rho - np.cos(np.pi*x) )
#plt.plot( x, rho )
plt.show()
