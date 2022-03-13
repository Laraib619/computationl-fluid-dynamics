import numpy as np
import matplotlib.pyplot as plt
import csv

## defining constants and Variables
n =int(input("enter the size of the grid:"))
l=1
count= 0
m=2
x =np.empty([n], dtype=float)
yn= np.empty([n], dtype=float)
Q =np.empty([n, 1], dtype=float)
Qn=np.empty([n, 1], dtype=float)

dx = l/(n-1)
eps =.00001
D = 2+(((m*l)**2)*(dx)**2)
print("d =",D)
# Generating Mesh
x[0] =0
for i in range (1,n):
    x[i] = x[i-1] +(dx)
#Array initialization
for i in range(1,n):
    Qn[i] = 0
    Q[i] = 0
#left boundary condition
Q[0] =1
Qn[0] =1
# Iteration loop starts #
err= 10
while err>eps:
    err=0
    for i in range(1,n-1):
        Qn[i] = (Q[i+1] + Q[i-1])/D
        yn = abs(Qn[i] - Q[i])
        err = max(yn,err)
        print("err" , err)

        Q[i] = Qn[i]
    count+= 1
    print('count = ',count)
    Qn[n-1] = 2*Q[n-2]/D
    Q[n-1] = Qn[n-1]

print("\n-------------------------------------------------------------\n")
print("\nThe values of Theta obtained after ", count, "number of iterations:\n")
# file writing
with open('1D heat fin.csv ', 'w') as f:
    fl = csv.writer(f)
    fl.writerow(["col1",  "i",  "kval"])
data = []
for i in range(0, n):
    print(" Theta ", i, " = ", Qn[i], "\n")
    with open('1D heat fin.csv a ', 'a') as f:
        fl = csv.writer(f)
        fl.writerow([" Theta ", i, Qn[i]])
    print("\n")
##  Plotting by matplotlib

fig, ax = plt.subplots(figsize=(10, 10))
ax.plot(Qn)
ax.set(title = "1D heat conduction fin",
       xlabel = "No. of grid points ",
       ylabel = "values of Theta")
plt.show()