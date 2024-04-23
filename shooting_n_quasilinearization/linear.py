import numpy as np
import matplotlib.pyplot as plt

def y(y_arr, x, h):
    idx = int(x / h)
    return y_arr[idx]

def FillMatrices(x_arr, y_arr, y0, y1, A, f, N, h):
    for k in range(1, N):
        xk = x_arr[k]
        yk = y(y_arr, xk, h)
        h2 = h**2

        A[k][k - 1] = 1.0 / h2

        A[k][k] = -2.0 / h2 - xk / 2.0 / np.sqrt(yk)

        A[k][k + 1] = 1.0 / h2 

        f[k] = xk / 2.0 * np.sqrt(yk)

    A[0][0] = 1.0    
    A[N][N] = 1.0

    f[0] = y0
    f[N] = y1
    return A, f

def linearMethod(x_arr, y_arr, y0, y1, A, f, N, h, eps):
    A, f = FillMatrices(x_arr, y_arr, y0, y1, A, f, N, h)
    new = np.linalg.solve(A, f)

    iteration = 1
    while np.max(np.abs(y_arr - new)) > eps:
        label = f"iteration = {iteration}"
        #plt.plot(x_arr, y_arr[:-1], label=label)
        y_arr = new
        A, f = FillMatrices(x_arr, y_arr, y0, y1, A, f, N, h)
        new = np.linalg.solve(A, f)
        iteration += 1

    return x_arr, y_arr, iteration

def main():
    y0 = 0
    y1 = 2
    eps = 1e-7
    h = 1e-4
    x_arr = np.arange(0.0, 1.0, h)
    N = x_arr.size
    y_arr = np.full(N + 1, 1.0)
    A = np.zeros((N + 1, N + 1))
    f = np.zeros(N + 1)

    x_arr, y_arr, iteration = linearMethod(x_arr, y_arr, y0, y1, A, f, N, h, eps)
    plt.plot(x_arr, y_arr[:-1])
    plt.title(f"Final linear solution on iteration = {iteration}")
    plt.grid()
    plt.xlabel("X")
    plt.xlabel("Y")
    plt.show()


if __name__ == "__main__":
  main()
