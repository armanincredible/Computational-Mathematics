import numpy as np
import matplotlib.pyplot as plt
import util

def euler_method(step_size, t_var, state_vector, differential_function):
    for j in range(1, np.shape(t_var)[0]):
        k1_vec = differential_function(t_var[j - 1], state_vector[j - 1])
        state_vector[j] = state_vector[j - 1] + step_size * k1_vec
    return t_var, state_vector

def rk4_method(step_size, t_var, state_vector, differential_function):
    for j in range(1, np.shape(t_var)[0]):
        k1_vec = differential_function(t_var[j - 1], state_vector[j - 1])
        k2_vec = differential_function(t_var[j - 1] + 0.5 * step_size, state_vector[j - 1] + 0.5 * step_size * k1_vec)
        k3_vec = differential_function(t_var[j - 1] + 0.5 * step_size, state_vector[j - 1] + 0.5 * step_size * k2_vec)
        k4_vec = differential_function(t_var[j - 1] + step_size, state_vector[j - 1] + step_size * k3_vec)
        state_vector[j] = state_vector[j - 1] + step_size * (1 / 6 * k1_vec + 2 / 6 * k2_vec + 2 / 6 * k3_vec + 1 / 6 * k4_vec)
    return t_var, state_vector

def adams_bashforth_method(step_size, t_var, state_vector, differential_function):
    n_steps = np.shape(state_vector)[0]
    f_vector = np.zeros((n_steps, np.shape(state_vector[0])[0]))
    f_vector[0] = differential_function(t_var[0], state_vector[0])
    state_vector[1] = state_vector[0] + step_size * f_vector[0]
    f_vector[1] = differential_function(t_var[1], state_vector[1])
    state_vector[2] = state_vector[1] + step_size * f_vector[1]
    f_vector[2] = differential_function(t_var[2], state_vector[2])
    for j in range(3, np.shape(t_var)[0]):
        state_vector[j] = state_vector[j - 1] + step_size * (23 / 12 * f_vector[j - 1] - 16 / 12 * f_vector[j - 2] + 5 / 12 * f_vector[j - 3])
        f_vector[j] = differential_function(t_var[j], state_vector[j])
    return t_var, state_vector

def get_differential_equation_function(A, B):
    return lambda t, state_vector: np.array(
        [A + state_vector[1] * state_vector[0] ** 2 - (B + 1) * state_vector[1], 
         B * state_vector[0] - state_vector[1] * state_vector[0] ** 2])

for b_value in range(2, 6):
    differential_equation = get_differential_equation_function(1, b_value)
    step_size = 0.001

    for method_name, method_func, color in [("Euler (1 Order)", euler_method, "r"),
                                     ("RK4 (4 Order)", rk4_method, "g"),
                                     ("Adams Bashforth(3 Order)", adams_bashforth_method, "b")]:
        
        time_steps, solution_vector = util.solve_differential_equation(step_size, 0, 100, np.array([1, 1]), 
                                                                        differential_equation, method_func)
        
        x1_values = np.array([solution_vector[i][0] for i in range(np.shape(solution_vector)[0])])
        x2_values = np.array([solution_vector[i][1] for i in range(np.shape(solution_vector)[0])])

        plt.plot(x1_values, x2_values, c=color, label=method_name)

    plt.legend()
    plt.xlabel('x1')
    plt.ylabel('x2')
    plt.title(f'B: {b_value} and Step Size: {step_size}')
    plt.show()
