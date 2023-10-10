import direct as direct
import iteration as iter
import scipy
import matplotlib.pyplot as plt
import numpy as np
import matrix as mx
import math

def plot(title, data):
    plt.figure(figsize=(15, 15))

    plt.xlabel("iteration number")
    plt.ylabel("log(abs_accurancy)")
    plt.yscale("log")

    plt.title(title)
    iterations = [i for i in range(len(data))]
    plt.plot(iterations, data, ".-")
    plt.show()

def plot_few(title, names, datas):
    plt.figure(figsize=(15, 15))

    plt.xlabel("iteration number")    
    plt.ylabel("log(abs_accurancy)")
    plt.yscale("log")

    plt.title(title)
    for ind in range(len(datas)):
        iterations = [i for i in range(len(datas[ind]))]
        plt.plot(iterations, datas[ind], ".-", label=names[ind])
        plt.legend()
    plt.show()


eps = 1e-6

mat = mx.GetMatrix()
b = mx.GetB()

x = direct.Gauss(mat, b)
assert(mx.norm(np.matmul(mat,x) - b) < eps)

x = direct.LUmethod(mat, b)
assert(mx.norm(np.matmul(mat,x) - b) < eps)

A = np.array([[4.0,-1.0,1.0],[2.0,5.0,2.0],[1.0,2.0,4.0]])
b = np.array([8.0,3.0,11.0])
x, data = iter.Jacobi(A, b)
plot("Jacobi", data)

x, data = iter.GaussSeidel(A, b)
plot("Gauss", data)

data_for_relaxation = []
names = []
for w in np.arange(1., 2.2, 0.2):
    w = math.ceil(w * 10) / 10
    x, data = iter.UpperRelaxation(A, b, w)
    data_for_relaxation.append(data)
    names.append("w=" + str(w))

plot_few("UpperRelaxation", names, data_for_relaxation)