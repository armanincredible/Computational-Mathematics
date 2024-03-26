import numpy as np


def solve_differential_equation(h_step, t_start, t_end, x0_vec, diff_function, method):
    assert t_start <= t_end
    assert h_step > 0

    n_steps = int((t_end - t_start) / h_step) + 1
    t = np.linspace(t_start, t_end, n_steps)

    x_vec = np.zeros((n_steps, np.shape(x0_vec)[0]))
    for i in range(np.shape(x0_vec)[0]):
        x_vec[0][i] = x0_vec[i]

    return method(h_step, t, x_vec, diff_function)
