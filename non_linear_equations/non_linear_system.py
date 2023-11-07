import math
import matplotlib.pyplot as plt
import numpy as np

def F(x, y):
    return [math.sin(x + 1) - y - 1.2, 2*x + math.cos(y) - 2]

def mpi_get_x(x, y):
    return (2 - math.cos(y)) / 2

def mpi_get_y(x, y):
    return math.sin(x + 1) - 1.2

def J(x, y):
    matrix = np.zeros((2, 2))
    matrix[0][0] = math.cos(x + 1)
    matrix[0][1] = -1
    matrix[1][0] = 2
    matrix[1][1] = -math.sin(y)
    return matrix

def res_plot(title, res, x, y):
    plt.grid()

    plt.title(title)
    plt.xlabel("X")
    plt.ylabel("Y")

    plt.plot(x, y, '.-', ms=8.0)
    plt.savefig(res)
    plt.clf()

def mpi_sycle_operation(x, y):
    return mpi_get_x(x, y), mpi_get_y(x, y)

def newton_sycle_operation(x, y):
    dF = J(x, y)
    x1 = x - (np.linalg.inv(dF) @ F(x, y))[0]
    y1 = y - (np.linalg.inv(dF) @ F(x, y))[1]
    return x1, y1

def simple_sycle(operation, start_values):
    x1 = start_values[0]
    y1 = start_values[1]
    x = x1 + 1
    y = y1 + 1
    xs = []
    ys = []
    while (abs(x1 - x) > 1e-3) and (abs(y1 - y) > 1e-3):
        x = x1
        y = y1
        xs.append(x)
        ys.append(y)
        x1, y1 = operation(x, y)
    print("res: ", x1, y1)
    return xs, ys

print("MPI:")

xs, ys = simple_sycle(mpi_sycle_operation, [5, 5])

res_plot("XY MPI", "img/system_mpi.png", xs, ys)

print("Newton:")

xs, ys = simple_sycle(newton_sycle_operation, [5, 5])
    
res_plot("XY Newton", "img/system_newton.png", xs, ys)

