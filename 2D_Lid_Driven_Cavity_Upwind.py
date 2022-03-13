import numpy as np
# import pandas as pd

import csv

if __name__ == '__main__':
    m=21
    i, j,k, Re = 400,400,400,400
    U = 1
    L = 1
    H = 1
    count = 0
    B= nu= dx= dy= dt= t = 0
    eps = 0.0001
    w, wnew, wdiff = np.empty([m,m] , dtype=float),np.empty([m,m] , dtype=float),np.empty([m,m] , dtype=float)
    Y, Ynew, Ydiff = np.empty([m,m] , dtype=float),np.empty([m,m] , dtype=float),np.empty([m,m] , dtype=float)
    u = np.empty([m,m] , dtype=float)
    v = np.empty([m,m] , dtype=float)
    eps = 0.0001
    dx = 1.0 / (m - 1)
    dy = 1.0 / (m - 1)
    dt = 0.001
    B = dx / dy
    step = t / dt
    cfl = U * dt / (dx)
    print("The Max cfl number is = ",cfl)
    for i in range(0, m):
        for j in range(0, m):
            Y[i][j] = 0
            Ynew[i][j] = 0
            w[i][j] = 0

    for j in range(1, m - 1):
        w[0][j] = -(2 * Y[1][j]) / (dx * dx)
    for i in range(0, m):
        w[i][0] = -(2 * Y[i][1]) / (dy * dy)
    for j in range(1, m - 1):
        w[m - 1][j] = -(2 * Y[m - 2][j]) / (dx * dx)
    for i in range(0, m):
        w[i][m - 1] = -(2 * Y[i][m - 2] + 2 * U * dy) / (dy * dy)

    e1 = 100
    while e1 > eps:
        # if e1 == 100:
        e1 = 0
        for i in range(1, m - 1):
            for j in range(1, m - 1):
                u[i][j] = (Y[i][j + 1] - Y[i][j - 1]) / (2 * dy)
                v[i][j] = -(Y[i + 1][j] - Y[i - 1][j]) / (2 * dx)

        #print("u =",i,j,u)
        #print("v =",i,j,v)
        for i in range (1,m-1):
            for j in range (1,m-1):
                if( u[i][j] >=0 and v[i][j]>= 0):
                    #print("x")
                    term1 = w[i][j] - dt * (
                            u[i][j] * (w[i][j] - w[i - 1][j]) / (dx) + v[i][j] * (w[i][j] - w[i][j - 1]) / (dy))

                if (u[i][j]<0 and v[i][j]>=0):
                    term1 = w[i][j] - dt * (
                                u[i][j] * (w[i + 1][j] - w[i][j]) / (dx) + v[i][j] * (w[i][j] - w[i][j - 1]) / (dy))

                if (u[i][j]>=0 and v[i][j]<0):
                    term1 = w[i][j] - dt * (
                                u[i][j] * (w[i][j] - w[i - 1][j]) / (dx) + v[i][j] * (w[i][j + 1] - w[i][j]) / (dy))
                if  u[i][j]<0 and v[i][j]<0:
                    term1 = w[i][j] - dt * (
                                u[i][j] * (w[i + 1][j] - w[i][j]) / (dx) + v[i][j] * (w[i][j + 1] - w[i][j]) / (dy))



                term2 = (dt / Re) * ((w[i + 1][j] - 2 * w[i][j] + w[i - 1][j]) / (dx * dx) + (
                                w[i][j + 1] - 2 * w[i][j] + w[i][j - 1]) / (dy * dy))

                wnew[i][j] = term1+ term2                    #term1 + term2

                wdiff[i][j] = abs(wnew[i][j] - w[i][j])

                e1 = max(wdiff[i][j], e1)

                w[i][j] = wnew[i][j]


        e2 = 100
        while e2> eps:
            e2=0
            # print("e2.", e2, eps)
            for i in range(1,m-1):
                for j in range (1,m-1):
                    Ynew[i][j] = (Y[i + 1][j] + Y[i - 1][j] + Y[i][j + 1] + Y[i][j - 1] + (dx * dx) * wnew[i][j]) / 4

                    Ydiff[i][j] = abs(Ynew[i][j] - Y[i][j])

                    e2 = max(Ydiff[i][j], e2)


                    Y[i][j] = Ynew[i][j]

        for j in range(1,m-1):
            w[0][j] = -(2 * Ynew[1][j]) / (dx * dx)

        for i in range (0,m):
            w[i][0] = -(2 * Ynew[i][1]) / (dy * dy)

        for j in range(1,m-1):
            w[m - 1][j] = -(2 * Ynew[m - 2][j]) / (dx * dx)

        for i in range(0,m):
            w[i][m - 1] = -(2 * Ynew[i][m - 2] + 2 * U * dy) / (dy * dy)

        t= t+ dt
        #print("count",count)
        count=count+1

    print("\n------------------------------------------------\n")
    print("The final values of u, v, w and Y after ",count," number of iterations are as follows: \n")

    for j in range(1,m-1):
        print("u",[j+1],"=",u[(m-1)//2][j])

    print("\n----------------------------------------------------\n")

    for i in range(1,m-1):
        print("v" ,[i+1],"=" ,v[i][(m-1)//2])

    f = open("2D_Lid_Driven_Cavity_Upwind_400Re_v_vel.txt", "a")
    f.write("Case: Re=100, 21*21 GP, dt=0.001, CFl= 0.128\n")

    for j in range(0, m):
        f.write(str(u[j]) + "\t" + str(u[(m - 1) // 2][j]) + "\n")

    f.write("--------------------------------------------------------------\n\n\n")

    for i in range(0, m):
        f.write(str(u[i]) + "\t" + str(v[i][(m - 1) // 2]) + "\n")

    print("File written succesfully")





