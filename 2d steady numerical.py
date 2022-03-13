


import numpy as np
import matplotlib.pyplot as plt
import csv
x= int(input("enter the no. of nodes in x-direction:"))
y= int(input("enter the no.of nodes in y- direction:"))
l=1
w=1
count = 0
Grid_x=[0]*x
Grid_y =[0]*y
dx = l/(x-1)
dy = w/(y-1)
eps=.0001
Q, Qnew1 ,Qdiff =  np.empty([x,y],dtype=float),np.empty([x,y],dtype=float),np.empty([x,y],dtype=float),


Grid_x[0]=0
Grid_y[0]=0
for i in range (1,x):
    Grid_x[i] = Grid_x[i-1]+dx
#print("Grid_x=",Grid_x)


for j in range(1,y):
    Grid_y[j] = Grid_y[j-1]+ dy

for i in range(0,x):
    for j in range(0,y):

        Q[i][j] = 0
        Qnew1[i][j]= 0
        Qdiff[i][j]=0

for i in range(0,x):
    Q[i][y-1] = 100
    Qnew1[i][y-1]=100

err=10000
while err>eps:
    err = 0
    for i in range(1,x-1):
        for j in range(1,y-1):
            Qnew1[i][j] = (Q[i+1][j]+Q[i-1][j]+Q[i][j+1]+Q[i][j-1])/4


            Qdiff[i][j] = abs(Qnew1[i][j]-Q[i][j])
            err = max(err,Qdiff[i][j])
            # print("err=",  err)
    for i in range(0,x):
        for j in range(0,y):

            Q[i][j]= Qnew1[i][j]

    count += 1
    print("count=",count)
print("\n-------------------------------------------------------------\n")
print("\nThe values of temp obtained after ", count, "number of iterations:\n")
# for i in range(0,x):
#     for j in range(0,y):
# print(Qnew1[5][5])


with open('file.csv 2d', 'w') as f:
    fl = csv.writer(f)
    fl.writerow(["col1", "i", "j", "kval"])
data = []

for i in range(0, x):
    dd = []
    for j in range(0, y):
        print(" Q ", i, j, " = ",Qnew1[i][j], "\n")
        dd.append(Qnew1[i][j])
        # f = open('data.txt 2dn', 'a')
        # f.write(f"{i} {j} {Qnew1[i][j]} \n")
        # f.close()
    data.append(dd)
print("\n")

import matplotlib.pyplot as plt
a = np.linspace(0, y-1, y)
b = a
X, Y = np.meshgrid(a, b)
Z = np.array(data)
fig, ax = plt.subplots(figsize=(6,6))
ax.contourf(Y,X,Z,cmap='RdGy')
ax.set(title = "2D steady heat conduction plate using GS",
       xlabel = "No. of grid points in x ",
       ylabel = "no. of grid points in y")
plt.show()







