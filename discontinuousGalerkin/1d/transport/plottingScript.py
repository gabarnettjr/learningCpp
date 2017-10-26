import numpy as np
import matplotlib.pyplot as plt
import time

u = float( np.loadtxt("u.txt") )
t = float( np.loadtxt("t.txt") )
dt = float( np.loadtxt("dt.txt") )
x = np.loadtxt("x.txt")
a = min(x)
b = max(x)
nTimesteps = int( np.loadtxt("nTimesteps.txt") )

plt.ion()
for i in range(nTimesteps+1) :
    rho = np.loadtxt('./snapshots/'+str(i).zfill(5)+'.txt')
    plt.plot( x, rho )
    plt.axis( [a,b,-1.2,1.2] )
    plt.title( '{0:02.3f}'.format(t) )
    plt.draw()
    plt.clf()
    t = t + dt

plt.ioff()
# plt.plot( x, rho - np.exp(-10*x**2) )
# plt.plot( x, rho - np.cos(np.pi*x) )
plt.plot( x, rho )
plt.show()
