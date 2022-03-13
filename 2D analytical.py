import numpy as np
import matplotlib.pyplot as plt
import math
import csv
x = int(input("enter the no of nodes in x-direction:"))
y = int(input("enter the number of nodes in y - direction:"))
n
L= 1
To = 0
W =1
count=0
Tc = 100
T1 = np.empty([x,y],dtype=float)
T2 = np.empty([x,y],dtype=float)
T3 = np.empty([x,y],dtype=float)
T4 = np.empty([x,y],dtype=float)
T5 = np.empty([x,y],dtype=float)
T = np.empty([x,y],dtype=float)
Tnew = np.empty([x,y],dtype=float)
Grid_x=[0]*x
Grid_y =[0]*y

dx = L/(x-1)
dy = L/(y-1)
# print("dx = ",dx)

Grid_x[0]=0
for i in range(1,x):
    Grid_x[i]= Grid_x[i-1]+dx

Grid_y[0]=0
for j in range(1,y):
    Grid_y[j]=Grid_y[j-1]+dy


##array initialisation

for i in range(0,x):
    for j in range(0,y):
        T1[i][j] =0
        T2[i][j] = 0
        T3[i][j] = 0
        T4[i][j] = 0
        T5[i][j] = 0
        T[i][j]  = 0
        Tnew[i][j] = 0

### boundary condition
for i in range(0,x):
    Tnew[i][0] = 0     #bottom wall

for i in range(0,x):
    Tnew[i][y-1] =100    ## Top wall

for j in range(1,y-1):
    Tnew[0][j] =  0      ## left wall

for j in range (1,y-1):
    Tnew[x-1][j] =0      ## right wall
# print("Grid_x",Grid_x)

i=0
n=1
while n<9:
    for i in range(1,x-1):
        for j in range(1,y-1):
            T1[i][j] = (1-(-1)**n)/(n)
            T2[i][j] = (math.sinh(n*math.pi*Grid_y[j]/L))
            T3[i][j] =(math.sin(n*math.pi*(Grid_x[i])/L))
            T4[i][j] =(math.sinh(n*math.pi*W/L))

            T5[i][j]= (T1[i][j]*T2[i][j]*T3[i][j]/T4[i][j])

            T[i][j]= (2* T5[i][j]*Tc)/(math.pi)

            Tnew[i][j] =(T[i][j]+Tnew[i][j])


    n=n+2
    count +=1
print("\n-------------------------------------------------------------\n")
print("\nThe values of temp obtained after ", count, "number of iterations:\n")

# print("T1 temp=",T1[5][5])
# print("T2 temp=",T2[5][5])
# print("T3 temp=",T3[5][5])
# print("T4 temp=",T4[5][5])
# print("T5 temp=",T5[5][5])


with open('file.csv 2d', 'w') as f:
    fl = csv.writer(f)
    fl.writerow(["col1", "i", "j", "kval"])
data = []

for i in range(0, x):
    dd = []
    for j in range(0, y):
        print(" T ", i, j, " = ",Tnew[i][j], "\n")
        dd.append(Tnew[i][j])
        f = open('2D analytical data.txt', 'a')
        f.write(f"{i} {j} {Tnew[i][j]} \n")
        f.close()
    data.append(dd)
print("\n")


a = np.linspace(0, y-1, y)
b = a
X, Y = np.meshgrid(a, b)
Z = np.array(data)
fig, ax = plt.subplots(figsize=(6,6))
ax.contourf(Y,X,Z,cmap='RdGy')
ax.set(title = "2D steady heat conduction plate analytical",
       xlabel = "No. of grid points in x ",
       ylabel = "no. of grid points in y")
plt.show()







