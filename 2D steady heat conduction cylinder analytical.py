import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import scipy.special as spl
m = int(input("enter the no of nodes in r-direction:"))
n = int(input("enter the number of nodes in z - direction:"))
k=0
R= 1
To = 0
Z=3
count=0
Tc = 100
L = np.empty([k],dtype=float)
T1 = np.empty([m,n],dtype=float)
T2 = np.empty([m,n],dtype=float)
T3 = np.empty([m,n],dtype=float)
T4 = np.empty([m,n],dtype=float)
T5 = np.empty([m,n],dtype=float)
T = np.empty([m,n],dtype=float)
Tnew = np.empty([m,n],dtype=float)
Grid_r= [0]*m
Grid_z = [0]*n

dr = R/(m-1)
dz = Z/(n-1)
# print("dx = ",dx)

Grid_r[0]=0
for i in range(1,m):
    Grid_r[i]= Grid_r[i-1]+dr

Grid_z[0]=0
for j in range(1,n):
    Grid_z[j]=Grid_z[j-1]+dz


##array initialisation

for i in range(0,m):
    for j in range(0,n):
        T1[i][j] =0
        T2[i][j] = 0
        T3[i][j] = 0
        T4[i][j] = 0
        T5[i][j] = 0
        T[i][j]  = 0
        Tnew[i][j] = 0

### boundary condition
#top wall
for i in range(0,m):
    Tnew[i][n-1] = 1

#bottom wall
for i in range(0,m):

    Tnew[i][0] = 0
#right wall
for j in range(0,n):
    Tnew[m-1][j] =0

#
L=np.loadtxt('roots.txt')

i=0

while k<7:
    for i in range(1,m-1):
        for j in range(1,n-1):
            #for k in range(1,6):

            T1[i][j] = 1/(L[k]*R)
            T2[i][j] = math.sinh((L[k]*Grid_z[j])*spl.j0(L[k]*Grid_r[i]))
            T3[i][j] = math.sinh((L[k]*Z)*spl.y0(L[k]*R))


            T4[i][j]=  ((T2[i][j])/(T3[i][j]))

            T[i][j]= (2* T1[i][j]*T4[i][j])

            Tnew[i][j] =(T[i][j]+Tnew[i][j])

    for j in range(1,m-1):

        Tnew[0][j]= Tnew[0][j]



    k=k+1
    count +=1
print("\n-------------------------------------------------------------\n")
print("\nThe values of temp obtained after ", count, "number of iterations:\n")

with open('file.csv 2d', 'w') as f:
    fl = csv.writer(f)
    fl.writerow(["col1", "i", "j", "kval"])
data = []

for i in range(0, m):
    dd = []
    for j in range(0, n):
        print(" T ", i, j, " = ",Tnew[i][j], "\n")
        dd.append(Tnew[i][j])
        f = open('2D analytical cyl data.txt', 'a')
        f.write(f"{i} {j} {Tnew[i][j]} \n")
        f.close()
    data.append(dd)
print("\n")


a = np.linspace(0, n-1, n)
b = a
X, Y = np.meshgrid(a, b)
Z = np.array(data)
fig, ax = plt.subplots(figsize=(6,6))
ax.contourf(Y,X,Z,cmap='RdGy')
ax.set(title = "2D steady heat conduction cylinder analytical",
       xlabel = "No. of grid points in R ",
       ylabel = "no. of grid points in Z")
plt.show()

