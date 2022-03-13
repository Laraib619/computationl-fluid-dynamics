
def main():

    m = int(input("Enter the size of square grid:"))
    print("m", m)
    i, j, k, Re = 100, 100, 100, 100
    U, L, H = 1, 1, 1
    count = 0  # U is lid velocity, L is y length, H is x length111
    B, nu, dx, dy, dt, t = 0, 0, 0, 0, 0, 0
    step, cfl, e1, e2, term1, term2, term3, t1, t2, eps = 0.0001, 0.0001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.001  # B is beta=1 nu is kinematic viscosity
    Grid_X = [0] * m
    Grid_Y = [0] * m
    w = [[0 for i in range(m)] for i in range(m)]
    wnew = [[0 for i in range(m)] for i in range(m)]
    wdiff = [[0 for i in range(m)] for i in range(m)]
    Y, Ynew, Ydiff = [[0 for i in range(m)] for i in range(m)], [[0 for i in range(m)] for i in range(m)], [
        [0 for i in range(m)] for i in range(m)]
    u, uf = [[0 for i in range(m)] for i in range(m)], [[0 for i in range(m)] for i in range(m)]
    v, vf = [[0 for i in range(m)] for i in range(m)], [[0 for i in range(m)] for i in range(m)]
    dx = 1.0 / (m - 1)
    dy = 1.0 / (m - 1)
    dt = 0.001
    B = dx / dy
    step = t / dt
    cfl = U * dt / (dx)

    print("The Reynolds number is = ", Re, "\n")
    print("The Max cfl number is = ", cfl, "\n")

    # //Grid generation:
    Grid_X[0] = 0
    for i in range(1, m):
        Grid_X[i] = Grid_X[i - 1] + dx

    Grid_Y[0] = 0
    for i in range(1, m):
        Grid_Y[i] = Grid_Y[i - 1] + dy

    # //Initialisation of arrays
    for i in range(0, m):

        for j in range(0, m):

            Y[i][j] = 0
            Ynew[i][j] = 0
            w[i][j] = 0


    # //Boundary Conditions

    for j in range(1, m - 1):
        w[0][j] = -(2 * Y[1][j]) / (dx * dx)  # //Left wall condition

    for i in range(0, m):
        w[i][0] = -(2 * Y[i][1]) / (dy * dy)  # //Bottom wall condition

    for j in range(0, m - 1):
        w[m - 1][j] = -(2 * Y[m - 2][j]) / (dx * dx)  # //Right wall condition

    for i in range(0, m):
        w[i][m - 1] = -(2 * Y[i][m - 2] + 2 * U * dy) / (dy * dy)  # //Top wall condition


    e1 = 100
    while e1>eps:
        e1=0

        for i in range(1, m - 1):

            for j in range(1, m - 1):
                u[i][j] = (Y[i][j + 1] - Y[i][j - 1]) / (2 * dy)
                v[i][j] = -(Y[i + 1][j] - Y[i - 1][j]) / (2 * dx)

        # //to calc face values of vel u & v by interpolation

        for i in range(1, m - 1):
            for j in range(1, m - 1):

                uf[i][j] = u[i - 1][j] + (u[i][j] - u[i - 1][j]) / 2

                vf[i][j] = v[i][j - 1] + (v[i][j] - v[i][j - 1]) / 2


        # //For points adjacent to boundary, CD is applied
        for j in range(1, m - 1):  # //All Points adjacent to left wall
            i = 1
            t1 = w[i][j] - dt * (
                    u[i][j] * (w[i + 1][j] - w[i - 1][j]) / (2 * dx) + v[i][j] * (w[i][j + 1] - w[i][j - 1]) / (
                    2 * dy))
            t2 = (dt / Re) * ((w[i + 1][j] - 2 * w[i][j] + w[i - 1][j]) / (dx * dx) + (
                    w[i][j + 1] - 2 * w[i][j] + w[i][j - 1]) / (dy * dy))
            wnew[i][j] = t1 + t2

            wdiff[i][j] = abs(wnew[i][j] - w[i][j])
            e1 = max(wdiff[i][j], e1)
            w[i][j] = wnew[i][j]

        for j in range(1, m - 1):  # //All Points adjacent to right wall

            i = m - 2
            t1 = w[i][j] - dt * (
                    u[i][j] * (w[i + 1][j] - w[i - 1][j]) / (2 * dx) + v[i][j] * (w[i][j + 1] - w[i][j - 1]) / (
                    2 * dy))
            t2 = (dt / Re) * ((w[i + 1][j] - 2 * w[i][j] + w[i - 1][j]) / (dx * dx) + (
                    w[i][j + 1] - 2 * w[i][j] + w[i][j - 1]) / (dy * dy))
            wnew[i][j] = t1 + t2

            wdiff[i][j] = abs(wnew[i][j] - w[i][j])
            e1 = max(wdiff[i][j], e1)
            w[i][j] = wnew[i][j]

        for i in range(1, m - 1):  # //All Points adjacent to bottom wall

            j = 1
            t1 = w[i][j] - dt * (
                    u[i][j] * (w[i + 1][j] - w[i - 1][j]) / (2 * dx) + v[i][j] * (w[i][j + 1] - w[i][j - 1]) / (
                    2 * dy))
            t2 = (dt / Re) * ((w[i + 1][j] - 2 * w[i][j] + w[i - 1][j]) / (dx * dx) + (
                    w[i][j + 1] - 2 * w[i][j] + w[i][j - 1]) / (dy * dy))
            wnew[i][j] = t1 + t2

            wdiff[i][j] = abs(wnew[i][j] - w[i][j])
            e1 = max(wdiff[i][j], e1)
            w[i][j] = wnew[i][j]

        for i in range(1, m - 1):  # //All Points adjacent to top wall

            j = m - 2
            t1 = w[i][j] - dt * (
                    u[i][j] * (w[i + 1][j] - w[i - 1][j]) / (2 * dx) + v[i][j] * (w[i][j + 1] - w[i][j - 1]) / (
                    2 * dy))
            t2 = (dt / Re) * ((w[i + 1][j] - 2 * w[i][j] + w[i - 1][j]) / (dx * dx) + (
                    w[i][j + 1] - 2 * w[i][j] + w[i][j - 1]) / (dy * dy))
            wnew[i][j] = t1 + t2

            wdiff[i][j] = abs(wnew[i][j] - w[i][j])
            e1 = max(wdiff[i][j], e1)
            w[i][j] = wnew[i][j]

            # //Quick Scheme for inner grip points
        for i in range(2, m - 2):

            for j in range(2, m - 2):

                if (uf[i][j] >= 0 and uf[i - 1][j] >= 0):
                    term1 = w[i][j] - (dt / dx) * u[i][j] * (
                            ((6 / 8) * w[i][j] + (3 / 8) * w[i + 1][j] - (1 / 8) * w[i - 1][j]) - (
                            (6 / 8) * w[i - 1][j] + (3 / 8) * w[i][j] - (1 / 8) * w[i - 2][j]))


                if (uf[i][j] < 0 and uf[i - 1][j] >= 0):

                    term1 = w[i][j] - (dt / dx) * u[i][j] * (
                            ((6 / 8) * w[i + 1][j] + (3 / 8) * w[i][j] - (1 / 8) * w[i + 2][j]) - (
                            (6 / 8) * w[i - 1][j] + (3 / 8) * w[i][j] - (1 / 8) * w[i - 2][j]))


                if (uf[i][j] >= 0 and uf[i - 1][j] < 0):

                    term1 = w[i][j] - (dt / dx) * u[i][j] * (
                            ((6 / 8) * w[i][j] + (3 / 8) * w[i + 1][j] - (1 / 8) * w[i - 1][j]) - (
                            (6 / 8) * w[i][j] + (3 / 8) * w[i - 1][j] - (1 / 8) * w[i + 1][j]))


                if (uf[i][j] < 0 and uf[i - 1][j] < 0):

                    term1 = w[i][j] - (dt / dx) * u[i][j] * (
                            ((6 / 8) * w[i + 1][j] + (3 / 8) * w[i][j] - (1 / 8) * w[i + 2][j]) - (
                            (6 / 8) * w[i][j] + (3 / 8) * w[i - 1][j] - (1 / 8) * w[i + 1][j]))


                if (vf[i][j] >= 0 and vf[i][j - 1] >= 0):

                    term2 = -(dt / dy) * v[i][j] * (
                            ((6 / 8) * w[i][j] + (3 / 8) * w[i][j + 1] - (1 / 8) * w[i][j - 1]) - (
                            (6 / 8) * w[i][j - 1] + (3 / 8) * w[i][j] - (1 / 8) * w[i][j - 2]))


                if (vf[i][j] < 0 and vf[i][j - 1] >= 0):

                    term2 = -(dt / dy) * v[i][j] * (
                            ((6 / 8) * w[i][j + 1] + (3 / 8) * w[i][j] - (1 / 8) * w[i][j + 2]) - (
                            (6 / 8) * w[i][j - 1] + (3 / 8) * w[i][j] - (1 / 8) * w[i][j - 2]))


                if (vf[i][j] >= 0 and vf[i][j - 1] < 0):

                    term2 = -(dt / dy) * v[i][j] * (
                            ((6 / 8) * w[i][j] + (3 / 8) * w[i][j + 1] - (1 / 8) * w[i][j - 1]) - (
                            (6 / 8) * w[i][j] + (3 / 8) * w[i][j - 1] - (1 / 8) * w[i][j + 1]))


                if (vf[i][j] < 0 and vf[i][j - 1] < 0):

                    term2 = -(dt / dy) * v[i][j] * (
                            ((6 / 8) * w[i][j + 1] + (3 / 8) * w[i][j] - (1 / 8) * w[i][j + 2]) - (
                            (6 / 8) * w[i][j] + (3 / 8) * w[i][j - 1] - (1 / 8) * w[i][j + 1]))


                term3 = (dt / Re) * ((w[i + 1][j] - 2 * w[i][j] + w[i - 1][j]) / (dx * dx) + (
                        w[i][j + 1] - 2 * w[i][j] + w[i][j - 1]) / (dy * dy))

                wnew[i][j] = term1 + term2 + term3

                wdiff[i][j] = abs(wnew[i][j] - w[i][j])
                e1 = max(wdiff[i][j], e1)

                w[i][j] = wnew[i][j]



        e2 = 100
        while e2>eps:
            e2=0

            for i in range(1, m - 1):

                for j in range(1, m - 1):

                    Ynew[i][j] = (Y[i + 1][j] + Y[i - 1][j] + Y[i][j + 1] + Y[i][j - 1] + (dx * dx) * wnew[i][
                        j]) / 4  # //Since B=1

                    Ydiff[i][j] = abs(Ynew[i][j] - Y[i][j])
                    e2 = max(Ydiff[i][j], e2)

                    Y[i][j] = Ynew[i][j]

            print("error2=",e2)





            # //Updating BC's for w

        for j in range(1, m - 1):
            w[0][j] = -(2 * Ynew[1][j]) / (dx * dx)  # //Left wall condition

        for i in range(0, m):
            w[i][0] = -(2 * Ynew[i][1]) / (dy * dy)  # //Bottom wall condition
        for j in range(1, m - 1):
            w[m - 1][j] = -(2 * Ynew[m - 2][j]) / (dx * dx)  # //Right wall condition
        for i in range(0, m):
            w[i][m - 1] = -(2 * Ynew[i][m - 2] + 2 * U * dy) / (dy * dy)  # //Top wall condition

        t = t + dt
        count += 1


    print("\n------------------------------------------------\n")
    print("The final values of u, v, Y after ", count, " number of iterations are as follows: \n")

    for j in range(0, m):
        print("u", j + 1, "=", u[(m - 1) // 2][j], "\n")

    print("\n----------------------------------------------------\n")
    for i in range(1, m - 1):

        print("v", i + 1, "=", v[i][(m - 1) // 2], "\n")


    # //Writing the results to a file

    f = open("2D_Lid_Driven_Cavity_Quick_100Re_v_vel.txt", "a")
    f.write("Case: Re=100, 129*129 GP, dt=0.001, CFl= 0.128\n")
    for i in range(0, m):

        for j in range(0, m):
            f.write(str(Grid_X[i]) + "\t" + str(Grid_Y[j]) + "\t" + str(Y[i][j]) + "\n")


    for j in range(0, m):
        f.write(str(Grid_Y[j]) + "\t" + str(u[(m - 1) // 2][j]) + "\n")


    f.write("--------------------------------------------------------------\n\n\n")


    for i in range(0, m):

        f.write(str(Grid_X[i]) + "\t" + str(v[i][(m - 1) // 2]) + "\n")

    print("File written successfully\n")
    f.close()

    import matplotlib.pyplot as plt
    import numpy
    a = numpy.linspace(0, m-1, m)
    b = a
    X, Y = numpy.meshgrid(a, b)
    Z = numpy.array(v)
    fig, ax = plt.subplots(figsize=(6,6))
    ax.contourf(X,Y,Z,cmap='RdGy')
    plt.show()



main()