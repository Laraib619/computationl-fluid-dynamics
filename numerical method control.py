import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc ('font',family='serif')
rc('lines',linewidth = 1.5)
rc('font',size = 16)
plt.rc('legend',**{'fontsize':14})

#######################################################################
#some definitions
#######################################################################
Nx = 401
Xmax= 10
Xmin = -10
dx = (Xmax- Xmin)/(Nx-1)
x = np.linspace(Xmin,Xmax,Nx)
a=1
lam=0.59
dt = .04
T =0
T_end = 5
Nt= int(T_end/dt)
t= np.linspace(0.,T_end,Nt+1)
   #c=.8
CFL = lam*(dt/dx)
U = np.zeros((Nt+1,Nx))
U[0,:] = max(1-abs(x))
Uex = U[0,:]
#############################################################################
#solve equaton using upwind scheme
#############################################################################
for n in range(0,Nt):
    if (lam>0.):
        for i in range(1,Nx):
            U[n+1,i] = U[n,i] - CFL*(U[n,i]-U[n,i-1])
        U[n+1,0] = U[n+1,Nx-1]
    else:
        for i in range(0,Nx-1):
            U[n+1,i]= U[n,i]- CFL*(U[n,i+1]-U[n,i])
        U[n+1,Nx-1]= U[n,0]
    #print(U)
##############################################################################
#compute exact solution
##############################################################################
    d =(lam*(n+1)*dt)
    Uex= (np.mod(x-d+Xmax,4)-Xmax)
    errL1= U-Uex
    errl2= np.linalg.norm(errL1)
##############################################################################
    if (n==0):  fig ,ax = plt.subplots(figsize=(5.5,4))
    plt.clf()
    plt.plot(x,U[n+1,:])
    plt.scatter(x,Uex,marker='o',facecolors='none',color= 'k',linewidth=1)
    plt.gca().legend(('upwind scheme (CFL='+str(CFL)+')','exact solution'))
    plt.axis([Xmin,Xmax,0,1.4])
    plt.title('t='+str(dt*(n+1)),fontsize = 16)
    plt.xlabel('x',fontsize=18)
    plt.ylabel('U',fontsize=18)
    plt.subplots_adjust(left= 0.2)
    plt.subplots_adjust(bottom = 0.18)
    plt.draw()
    plt.pause(0.001)

plt.show()
