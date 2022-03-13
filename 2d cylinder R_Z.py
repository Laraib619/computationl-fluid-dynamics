import numpy as np
import matplotlib.pyplot as plt
import csv

m = int(input("Enter the no. of nodes in R-direction:"))
n = int(input("Enter the no. of nodes in Z-direction:"))
r= 1
z =3
count= 0
Grid_r= [0]*m
Grid_z = [0]*n
eps = .0001
Q= np.empty([m,n],dtype=float)
Qn = np.empty([m,n],dtype=float)
Qdiff = np.empty([m,n],dtype=float)
dr= float(r/(m-1))
dz= float(z/(m-1))
print('dr=',dr)
print('dz=',dz)

Grid_r[0]= 0
Grid_z[0]= 0
for i in range(1,m):
    Grid_r[i] = Grid_r[i-1]+ dr

for j in range(1,n):
    Grid_z[i] = Grid_z[i-1] +dz

for i in range(0,m):
    for j in range(0,n):
        Q[i][j] = 0
        Qn[i][j] = 0
        Qdiff[i][j]=0
#top wall
for i in range(0,m):
    Q[i][n-1]= 1
    Qn[i][n-1] = 1

#bottom wall
for i in range(0,m):
    Q[i][0] = 0
    Qn[i][0] = 0
#right wall
for j in range(0,n):
    Q[m-1][j] =0
    Qn[m-1][j] =0



err =1000
while err>eps:
    err=0
    for i in range(1,m-1):
        for j in range(1,n-1):

            T1 = 1/(2*(1+(dr*dr)/(dz*dz)))

            T2 = Q[i+1][j]+Q[i-1][j]+(dr*dr)/(dz*dz)*(Q[i][j+1]+Q[i][j-1])+dr/(2*Grid_r[i])*(Q[i+1][j]-Q[i-1][j])
            Qn[i][j] = (T1*T2)



            Qdiff[i][j] = abs(Qn[i][j]-Q[i][j])
            err = max(err,Qdiff[i][j])
            #Q[i][j] = Qn[i][j]

    for i in range(0,m):
        for j in range(0,n):
            Q[i][j] = Qn[i][j]

    for j in range(1,n-1):
        Q[0][j]=Q[1][j]
        Qn[0][j]= Qn[0][j]
    count+=1
    print(count)

print("\n-------------------------------------------------------------\n")
print("\nThe values of temp obtained after ", count, "number of iterations:\n")
with open('file.csv 2d cylinder', 'w') as f:
    fl = csv.writer(f)
    fl.writerow(["col1", "i", "j", "kval"])
data = []

for i in range(0, m):
    dd = []
    for j in range(0, n):
        print(" Q ", i, j, " = ",Qn[i][j], "\n")
        dd.append(Qn[i][j])
        f = open('data.txt 2dn cylinder', 'a')
        f.write(f"{i} {j} {Qn[i][j]} \n")
        f.close()
    data.append(dd)
print("\n")

import matplotlib.pyplot as plt
a = np.linspace(0, n-1, n)
b = a
X, Y = np.meshgrid(a, b)
Z = np.array(data)
fig, ax = plt.subplots(figsize=(6,6))
ax.contourf(Y,X,Z,cmap='RdGy')
ax.set(title = "2D steady heat conduction cylinder using GS",
       xlabel = "No. of grid points in R ",
       ylabel = "no. of grid points in Z")
plt.show()

