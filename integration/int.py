import numpy as np
import pandas as pd
from scipy import integrate as sp_integrate
import matplotlib.pyplot as plt
from IPython.display import display, Markdown

def get_constant_step(x):
    step_sizes = np.diff(x)
    assert np.allclose(step_sizes, step_sizes[0]), "x step must be constant"
    return step_sizes[0]

def integrate_trapezoid(x, y):
    h = get_constant_step(x)
    return h * (y[0]/2 + y[-1]/2 + np.sum(y[1:-1]))

def integrate_richardson(x, y, int_function, order):
    val_2h = int_function(x[::2], y[::2])
    val_h = int_function(x, y)
    return val_h + (val_h - val_2h) / (2**order - 1)

def integrate_simpson(x, y):
    h = get_constant_step(x)
    weights = np.full(len(x), 4)
    weights[::2] = 2
    weights[0] = weights[-1] = 1
    return 1/3 * h * np.sum(weights * y)


"""-----------------------------------------------------SPLINES-----------------------------------------------------"""

def solve_tridiagonal(matrix, vector):
    shape = matrix.shape
    dimension = len(vector)
    assert shape[0] == dimension and shape[1] == dimension, "Dimensions must be equal"

    lower_diag = np.diag(matrix, -1)
    main_diag = np.diag(matrix, 0).copy()
    upper_diag = np.diag(matrix, 1)
    right_side = vector.copy()

    for i in range(1, dimension):
        w = lower_diag[i-1] / main_diag[i-1]
        main_diag[i] -= w * upper_diag[i-1]
        right_side[i] -= w * right_side[i-1]

    solution = np.zeros(dimension)
    solution[dimension-1] = right_side[dimension-1] / main_diag[dimension-1]

    for i in reversed(range(0, dimension-1)):
        solution[i] = (right_side[i] - upper_diag[i] * solution[i+1]) / main_diag[i]

    difference = np.matmul(matrix, solution) - vector
    assert np.linalg.norm(difference, ord=np.inf) < 1e-9, 'The solution is incorrect'

    return solution

def cubic_spline_interpolation(x_values, y_values):
    assert len(x_values) == len(y_values), "x and y must have the same size"
    dimension = len(x_values)

    h_values = np.diff(x_values)

    mu_values = h_values[:-1] / (h_values[:-1] + h_values[1:])

    lambda_values = 1 - mu_values
    lambda_values[0] = 0

    order = len(x_values)
    divided_diff = np.zeros((order, dimension))
    divided_diff[0] = y_values

    for i in range(1, order):
        prev_values = divided_diff[i-1]
        divided_diff[i, :order-i] = [(prev_values[m+1] - prev_values[m]) / (x_values[m+i] - x_values[m])
                                      for m in range(0, order - i)]

    second_order_diff = divided_diff[2]

    d_values = np.zeros(dimension)
    d_values[1:-1] = 6 * second_order_diff[:dimension-2]

    # Construct lambda-mu matrix
    spline_matrix = np.zeros((dimension, dimension))
    np.fill_diagonal(spline_matrix, 2)
    np.fill_diagonal(spline_matrix[1:], mu_values)
    np.fill_diagonal(spline_matrix[:, 1:], lambda_values)

    m_values = solve_tridiagonal(spline_matrix, d_values)

    def cubic_spline(arg):
        m = np.searchsorted(x_values, arg)
        m = np.clip(m, 1, dimension-1)
        h = x_values[m] - x_values[m-1]
        t = (arg - x_values[m-1]) / h
        t2 = t * t
        t3 = t2 * t
        a = m_values[m-1] * (x_values[m] - arg)**3 / (6 * h)
        b = m_values[m] * (arg - x_values[m-1])**3 / (6 * h)
        c = (y_values[m-1] - m_values[m-1] * h**2 / 6) * (x_values[m] - arg) / h
        d = (y_values[m] - m_values[m] * h**2 / 6) * (arg - x_values[m-1]) / h
        val = a + b + c + d
        return val

    return lambda arg: cubic_spline(arg)



def integrate_n_print(x, y):
    plt.scatter(x, y)
    plt.grid()
    plt.show()

    trapz = integrate_trapezoid(x, y)
    ref_trapz = np.trapz(x=x, y=y)
    richardson = integrate_richardson(x, y, integrate_trapezoid, 2)
    simpson = integrate_simpson(x, y)
    ref_simpson = sp_integrate.simpson(x=x, y=y)

    print(f"Richardson     : {richardson}")
    print(f"Simpson        : {simpson}")
    print(f"Simpson (ref)  : {ref_simpson}")
    print(f"Trapezoid      : {trapz}")
    print(f"Trapezoid (ref): {ref_trapz}")

def main():
    y = np.array([1.000000, 0.989616, 0.958851, 0.908852, 0.841471, 0.759188, 0.664997, 0.562278, 0.454649])
    x = np.arange(0, 2 + 0.25, 0.25)

    integrate_n_print(x, y)

    print("Now splines:")

    k = 100
    x_1 = np.array([-0.4, -0.1, 0.2, 0.5, 0.8])
    y_1 = np.array([1.9823, 1.6710, 1.3694, 1.0472, 0.64360])

    interp = cubic_spline_interpolation(x_1, y_1)

    lin = np.linspace(min(x_1), max(x_1), 501)
    plt.plot(lin, interp(lin))
    plt.scatter(x_1, y_1)
    plt.grid()

    equ = lambda x_1: np.sin(k * x_1) * interp(x_1)
    x_1 = lin
    y_1 = equ(lin)
    integrate_n_print(x_1, y_1)

if __name__ == "__main__":
    main()
