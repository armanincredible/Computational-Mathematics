import numpy as np
import matplotlib.pyplot as plt

# Solve task XI.9.3(a)

def d2y(x, y):
    return x * np.sqrt(y)

def f(x, f0, df, h):
    f_new = f0 + df(x, f0) * h

    eps = 1e-4
    for i in range(1000):
        if abs(f_new - f0) <= eps:
            return f_new
        f0 = f_new
        f_new = f0 + df(x, f0) * h

    print("Too many iterations!")
    return f_new

def calculateY(x_arr, y0, alpha, h):
    y = y0
    dy = alpha

    ys = []

    for x in x_arr:
        ys.append(y)
        y += dy * h
        dy = f(x, dy, d2y, h)

    #plt.plot(x_arr, ys)

    return ys

def shootingMethod(x_arr, y0, y1, calc_y, alpha, h_alpha, eps, h):
    y_new = calc_y(x_arr, y0, alpha, h)
    F = y_new[-1] - y1
    while abs(F) > eps:
        y_new = calc_y(x_arr, y0, alpha + h_alpha, h)
        dF = (y_new[-1] - y1 - F) / h_alpha

        alpha = alpha - F/dF

        y_new = calc_y(x_arr, y0, alpha, h)
        F = y_new[-1] - y1
    return y_new, alpha

def main():
    h = 4e-5
    x_arr = np.arange(0.0, 1.0, h)
    y0 = 0
    y1 = 2
    alpha = 0
    h_alpha = 1e-2
    eps = 1e-7

    y_new, alpha = shootingMethod(x_arr, y0, y1, calculateY, alpha, h_alpha, eps, h)

    plt.grid()
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.plot(x_arr, y_new, label="Solution")
    plt.title(f"Final solution of Differential Equation, alpha = {alpha:.4f}")
    plt.show()

if __name__ == "__main__":
  main()