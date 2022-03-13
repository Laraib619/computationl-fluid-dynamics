# -*- coding: utf-8 -*-
"""
Created on Thu May 13 05:04:16 2021

@author: Laraib Quamar
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
if __name__ == '__main__':

    i = j = x = y = t = 20
    count = 0
    alpha = 0.001
    C1 = C2 = lx = 1
    ly = 1
    dt = 10
    eps = 0.0001
                  #alpha is diffusivity
    print("Enter the no of grid point in x direction: \n")
    x = input()
    print("Enter the no of grid point in y direction: \n")
    y = input()
    x = int(x)
    y = int(x)
    A,M, AX, B ,K = np.empty([x,1],dtype=float),np.empty([x,1],dtype=float),np.empty([x,1],dtype=float),np.empty([x,1],dtype=float),np.empty([x,1],dtype=float)
    C , CL , D , L = np.empty([y,1],dtype=float),np.empty([y,1],dtype=float),np.empty([y,1],dtype=float),np.empty([y,1],dtype=float)

    T , Tnew1 , Tnew2 , diff =  np.empty([x,y],dtype=float),np.empty([x,y],dtype=float),np.empty([x,y],dtype=float),np.empty([x,y],dtype=float)

    dx = lx / (x - 1)
    dy = ly / (y - 1)
    step = t / dt
    C1 = alpha * dt / (dx * dx)  # 1st CFL number
    C2 = alpha * dt / (dy * dy)  # 2nd CFL number
    print("\nThe CFL number is ", C1, " and ", C2)

    # //Initialising of matrices
    for i in range(0, x):
        for j in range(0, y):
            T[i][j] = 300
            Tnew1[i][j] = 300
            Tnew2[i][j] = 300

    # Setting up BC's ==> Here, Drichlet BC is used T(dependant variable) at all walls known.

    for j in range(0, y):  # Left wall
        Tnew1[0][j] = 400
        Tnew2[0][j] = 400
        T[0][j] = 400

    for i in range(1, x - 1):             #Top Wall
        Tnew1[i][y - 1] = 600
        Tnew2[i][y - 1] = 600
        T[i][y - 1] = 600

    for j in range(0, y):  # Right wall

        Tnew1[x - 1][j] = 400
        Tnew2[x - 1][j] = 400
        T[x - 1][j] = 400

    for i in range(1, x - 1):  # Bottom wall

        Tnew1[i][0] = 600
        Tnew2[i][0] = 600
        T[i][0] = 600

    A[0] = 0
    A[x - 1] = 0
    for i in range(1, x - 1):
        A[i] = (alpha * dt) / (2 * dx * dx)

    M[0] = 0
    M[x - 1] = 0

    for i in range(1, x - 1):
        M[i] = (alpha * dt) / (2 * dx * dx)

    C[0] = 0
    C[y - 1] = 0
    for j in range(1, y - 1):
        C[j] = (alpha * dt) / (2 * dy * dy)

    CL[0] = 0
    CL[y - 1] = 0

    for j in range(1, y - 1):
        CL[j] = (alpha * dt) / (2 * dy * dy)

    # Loop to perform iteration
    err: int=0
    # count = 0
    print(A,"AUUUU")
    err = 1000000000
    while err > eps:
        err = 0
        print("errr",err,"epsss",eps)

        # 1st Sweeping
        for j in range(1, y - 1):

            B[0] = 1
            B[x - 1] = 1
            XX = []
            for i in range(1, x - 1):
                B[i] = -(1 + alpha * dt / (dx * dx))
                # XX = B

            # XX = B
            K[0] = T[0][j]
            K[x - 1] = T[x - 1][j]

            for i in range(1, x - 1):
                K[i] = -(T[i][j] + (alpha * dt / (2 * dy * dy)) * (T[i][j + 1] - 2 * T[i][j] + T[i][j - 1]))

            # //Thomas Algo
            # B = XX

            for i in range(1, x-1):

                B[i] = B[i] - ((M[i] * A[i - 1]) / B[i - 1])

            for i in range(1, x):
                K[i] = K[i] - (K[i - 1] * M[i] / B[i - 1])

            Tnew1[x - 1][j] = K[x - 1] / B[x - 1]
            for i in range(x - 2, 0, -1):
                Tnew1[i][j] = (K[i] - A[i] * Tnew1[i + 1][j]) / B[i]

        # //2nd Sweeping

        for i in range(1, x - 1):
            D[0] = 1
            D[y - 1] = 1
            for j in range(1, y - 1):

                D[j] = -(1 + alpha * dt / (dy * dy))

            L[0] = T[i][0]
            L[y - 1] = T[i][y - 1]
            for j in range(1, y - 1):
                L[j] = -(Tnew1[i][j] + (alpha * dt / (2 * dx * dx)) * (
                            Tnew1[i + 1][j] - 2 * Tnew1[i][j] + Tnew1[i - 1][j]))

            # //Thomas Algo
            for j in range(1, y - 1):

                D[j] = D[j] - (CL[j] * C[j - 1] / D[j - 1])

            for j in range(1, y):
                L[j] = L[j] - (L[j - 1] * CL[j] / D[j - 1])

            Tnew2[i][y - 1] = L[y - 1] / D[y - 1]
            for j in range(y - 2, 0, -1):
                Tnew2[i][j] = (L[j] - C[j] * Tnew2[i][j + 1]) / D[j]

        for i in range(1, x - 1):
            for j in range(1, y - 1):
                diff[i][j] = abs(Tnew2[i][j] - T[i][j])
                print("diff::::",diff[i][j])
                err = max(diff[i][j], err)

        print("err", err)

        count += 1

        for i in range(0, x):

            for j in range(0, y):
                T[i][j] = Tnew2[i][j]



    print("\n-------------------------------------------------------------\n")
    print("\nThe values of temp obtained after ", count, "number of iterations:\n")

    with open('file.csv', 'w') as f:
        fl = csv.writer(f)
        fl.writerow(["col1","i","j","kval"])
    data = []
    for i in range(0, x):
        dd = []
        for j in range(0, y):
            print(" Tnew(2) ", i, j, " = ", Tnew2[i][j], "\n")
            dd.append(Tnew2[i][j])
            f = open('data.txt', 'a')
            f.write(f"{i} {j} {Tnew2[i][j]} \n")
            f.close()
        data.append(dd)
    print("\n")


    a = np.linspace(0, x - 1, x)
    b = a
    X, Y = np.meshgrid(a, b)
    Z = np.array(data)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.contourf(X, Y, Z, cmap='RdGy')
    plt.show()