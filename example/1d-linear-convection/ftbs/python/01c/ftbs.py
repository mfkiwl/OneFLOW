# Remember: comments in python are denoted by the pound sign
import numpy                       #here we load numpy
from matplotlib import pyplot      #here we load matplotlib
import time, sys                   #and load some utilities


def set_fuction(u, ct) :
    for i in range(nx):
        xm = i * dx
        if 0.5 <= xm - ct <= 1:
            u[i] = 2
        else:
            u[i] = 1

nx = 41  # try changing this number from 41 to 81 and Run All ... what happens?
dx = 2 / (nx-1)
nt = 25    #nt is the number of timesteps we want to calculate
dt = .025  #dt is the amount of time each timestep covers (delta t)
c = 1      #assume wavespeed of c = 1

total_t = dt * nt
utheory = numpy.ones(nx)
x = numpy.linspace(0, 2, nx)
print('x=',x)
u = numpy.ones(nx)      #numpy function ones()

set_fuction(u,0)
set_fuction(utheory,c * total_t)
        
print('int(.5 / dx)=',int(.5 / dx))
print('int(1 / dx + 1)=',int(1 / dx + 1))
print(u)
print(utheory)

un = numpy.ones(nx) #initialize a temporary array

for n in range(nt):  #loop for values of n from 0 to nt, so it will run nt times
    un = u.copy() ##copy the existing values of u into un
    for i in range(1, nx): ## you can try commenting this line and...
    #for i in range(nx): ## ... uncommenting this line and see what happens!
        u[i] = un[i] - c * dt / dx * (un[i] - un[i-1])
        
pyplot.plot(x, utheory, "k-", linewidth=1.0, label="Exact solution")
pyplot.scatter(x, u, facecolor="none", edgecolor="blue", s=20, linewidths=0.5, label="FTBS scheme")
pyplot.xlabel("$x$")
pyplot.ylabel("$u$")
pyplot.title("Solution field")
pyplot.legend()
pyplot.tight_layout()
pyplot.show();