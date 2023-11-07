import math
import matplotlib.pyplot as plt

def derivative(func, x, h):
    return (func(x + h) - func(x)) / h

def function(x):
    return 2.0 * math.log10(x) - x / 2.0 + 1

def newton_operation(x):
    return x - function(x)/derivative(function, x, 1e-4)

def mpi_opearion(x):
    return 4.0 * math.log10(x) + 2.0

def res_plot(res_file, title):
    plt.grid()

    plt.title(title)
    plt.ylabel("X")
    plt.xlabel("Step")

    plt.plot(xs, '.-', ms=5.0)
    plt.savefig(res_file)
    plt.clf()

def simple_sycle(operation, start_value):
    xs = []
    x1 = start_value
    x = start_value + 1 #usefull value to start iteration
    while abs(x1 - x) > 1e-6:
        print(x1)
        x = x1
        xs.append(x)
        x1 = operation(x)
    print("res: ", x1)
    return xs

print("MPI:")

xs = simple_sycle(mpi_opearion, 10)

res_plot("img/equation_mpi.png", "MPI method")

print("Newton:")

xs = simple_sycle(newton_operation, 10)

res_plot("img/equation_newton.png", "Newthon method")
