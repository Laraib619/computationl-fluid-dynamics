import numpy as np
import matplotlib.pyplot as plt
import csv


m=31





p = m*m
print("the total no.of grid points =", p)
i=j=k=Re= 40
U,L,H = 1,1,1
count=0
eps =.01
w,wnew,wdiff = np.empty([m,m] , dtype=float),np.empty([m,m] , dtype=float),np.empty([m,m] , dtype=float)
Y,Ynew,Ydiff = np.empty([m,m] , dtype=float),np.empty([m,m] , dtype=float),np.empty([m,m] , dtype=float)
u , v = np.empty([m,m] , dtype=float), np.empty([m,m] , dtype=float)
dx=1./(m-1)
dy = 1./(m-1)
dt=0.001
B=dx/dy
t=0
step=t/dt
cfl=U*dt/(dx)
print("The max cfl=" ,cfl)
for i in range(0,m):
    for j in range(0,m):
        Y[i][j] = 0
        Ynew[i][j] = 0
        w[i][j] = 0
for j in range(1,m-1):
    w[0][j] = -(2 * Y[1][j]) / (dx * dx)
for i in range(0,m):
    w[i][0] = -(2 * Y[i][1]) / (dy * dy)
for j in range(1,m-1):
    w[m - 1][j] = -(2 * Y[m - 2][j]) / (dx * dx)
for i in range(0,m):
    w[i][m - 1] = -(2 * Y[i][m - 2] + 2 * U * dy) / (dy * dy)
#loop to perform iteration
e1=100
while e1>eps:
    e1=0


    for i in range(1,m-1):
        for j in range(1,m-1):
            u[i][j] = (Y[i][j + 1] - Y[i][j - 1]) / (2 * dy)
            v[i][j] = -(Y[i + 1][j] - Y[i - 1][j]) / (2 * dx)
    for i in range(1,m-1):
        for j in range(1,m-1):
            term1 = w[i][j] - dt * (
                        u[i][j] * (w[i + 1][j] - w[i - 1][j])/(2*dx)+v[i][j]*(w[i][j + 1]- w[i][j - 1]) / (
                            2 * dy))
            term2 = (dt / Re) * ((w[i + 1][j] - 2 * w[i][j] + w[i - 1][j]) / (dx * dx) + (
                        w[i][j + 1] - 2 * w[i][j] + w[i][j - 1]) / (dy * dy))
            wnew[i][j] = term1 + term2

            wdiff[i][j] = abs(wnew[i][j] - w[i][j])
            e1 = max(wdiff[i][j], e1)

            w[i][j] = wnew[i][j]


        print("e1 =",e1)
    e2=100
    while e2>eps:
        e2=0

        for i in range(1,m-1):
            for j in range(1,m-1):
                Ynew[i][j] = (Y[i + 1][j] + Y[i - 1][j] + Y[i][j + 1] + Y[i][j - 1] + (dx * dx) * wnew[i][
                    j]) / 4
                B = 1

                Ydiff[i][j] = abs(Ynew[i][j] - Y[i][j])
                e2 = max(Ydiff[i][j], e2)

                Y[i][j] = Ynew[i][j]
                #if e2 == eps:
            count += 1  #break
            print("error 2",e2,"no.of iter",count)


    for j in range(1,m-1):
        w[0][j] = -(2 * Ynew[1][j]) / (dx * dx)
    for i in range(0,m):
        w[i][0] = -(2 * Ynew[i][1]) / (dy * dy)
    for j in range(1,m-1):
        w[m - 1][j] = -(2 * Ynew[m - 2][j]) / (dx * dx)
    for i in range(0,m):
        w[m - 1][j] = -(2 * Ynew[m - 2][j]) / (dx * dx)
    t= t+dt
    count +=1


print("\n------------------------------------------------\n")
print("\nThe final values of u, v, w and Y after ", count, "number of iterations:\n")
for i in range(1,m-1):
    print("v",i+1,"=",v[i][15])
print("\n----------------------------------------------------\n")
for j in range(1,m-1):
    print("u", j + 1, "=", u[15][j])


f = open("2d_lid_cd.txt", "a")
f.write("Case:m=31, 31*31 GP, dt=0.001, CFl= 0.128\n")
for j in range(1, m-1):
    f.write("u" + str(j) + " = " + str(u[15][j]) + "\n")

f.write("--------------------------------------------------------------\n\n\n")

for i in range(1, m-1):
    f.write("v" + str(i) + " = " + str(v[i][15]) + "\n")

f.write("--------------------------------------------------------------\n\n\n")

print("File written successfully\n")
f.close()

















