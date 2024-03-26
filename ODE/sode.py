import numpy as np
import matplotlib.pyplot as plt
import util

def solve_non_linear_eq(get_jacobian, get_sys_equation, initial_vector, stop_iter_epsilon):
    n_iters = 0
    k_vector_next = np.copy(initial_vector)
    k_vector = np.copy(initial_vector) + np.ones(len(initial_vector)) * 2 * stop_iter_epsilon

    while (np.abs(np.linalg.norm(k_vector, ord=np.inf) - np.linalg.norm(k_vector_next, ord=np.inf)) > stop_iter_epsilon
           and n_iters <= 30000):
        k_vector = np.copy(k_vector_next)
        k_vector_next = k_vector - np.dot(np.linalg.inv(get_jacobian(k_vector)), get_sys_equation(k_vector))
        n_iters += 1

    return k_vector_next, n_iters

def system_equation(t, state_vector):
    return np.array([
        77.27 * (state_vector[1] + state_vector[0] * (1 - 8.375 * 10 ** (-6) * state_vector[0] - state_vector[1])),
        1.0 / 77.27 * (state_vector[2] - (1 + state_vector[0]) * state_vector[1]),
        0.161 * (state_vector[0] - state_vector[2]),
    ])

def jacobian(state_vector):
    return np.array([
        [77.27 * (1 - state_vector[0] * 2 * 8.375 * 10 ** (-6) - state_vector[1]), 77.27 * (1 - state_vector[0]), 0],
        [-1 / 77.27 * state_vector[1], 1 / 77.27 * (-1 - state_vector[0]), 1 / 77.27],
        [0.161, 0, -0.161],
    ])

def gear_jacobian_func(step_size, get_jacobian_func, state_vector):
    return lambda x: np.eye(np.shape(state_vector)[0]) - 6 * step_size / 11 * get_jacobian_func(x)

def gear_sys_equation_func(step_size, function, state_vector, t, x_free_vec):
    return lambda x: x - 6 * step_size / 11 * system_equation(t, x) + x_free_vec

def gear_method(step_size, t, state_vector, diff_function):
    n_steps = np.shape(state_vector)[0]
    f_vec = np.zeros((n_steps, np.shape(state_vector[0])[0]))

    f_vec[0] = diff_function(t[0], state_vector[0])
    state_vector[1] = state_vector[0] + step_size * f_vec[0]

    f_vec[1] = diff_function(t[1], state_vector[1])
    state_vector[2] = state_vector[1] + step_size * f_vec[1]

    f_vec[2] = diff_function(t[2], state_vector[2])

    for j in range(3, np.shape(t)[0]):
        x_free_vec = -18 / 11 * state_vector[j - 1] + 9 / 11 * state_vector[j - 2] - 2 / 11 * state_vector[j - 3]
        x_start_vec = state_vector[j - 1] + step_size * (23 / 12 * f_vec[j - 1] - 16 / 12 * f_vec[j - 2] + 5 / 12 * f_vec[j - 2])
        state_vector[j], _ = solve_non_linear_eq(
            gear_jacobian_func(step_size, jacobian, state_vector[j]),
            gear_sys_equation_func(step_size, system_equation, state_vector[j], t[j], x_free_vec),
            x_start_vec, 1e-6
        )

        f_vec[j] = diff_function(t[j], state_vector[j])

    return t, state_vector

t_step, solution_vec = util.solve_differential_equation(0.01, 0, 800, np.array([4, 1.1, 4]), system_equation, gear_method)

figure = plt.figure(figsize=(16, 7), facecolor='#F5F5F5')
axes = figure.add_subplot(1, 1, 1)

axes.xaxis.set_major_locator(plt.MaxNLocator(10))
axes.minorticks_on()
axes.grid(which='major', linewidth=2, color='#919191')
axes.grid(which='minor', linestyle=':')

x1_values = np.array([solution_vec[i][0] for i in range(np.shape(solution_vec)[0])])
x2_values = np.array([solution_vec[i][1] for i in range(np.shape(solution_vec)[0])])
x3_values = np.array([solution_vec[i][2] for i in range(np.shape(solution_vec)[0])])

axes.plot(t_step, x1_values, c='g', label='x1')
axes.plot(t_step, x2_values, c='r', label='x2')
axes.plot(t_step, x3_values, c='b', label='x3')

axes.legend()
axes.set_xlabel('t')
axes.set_ylabel('value')

plt.show()
