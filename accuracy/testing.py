import math
import matplotlib.pyplot as plt
import seaborn as sns

import accuracy as acc
import functions as fcs

MaxPowerOfDelta = 21
ConstantX = 10

ApproxDerivativeFunctions = [acc.approx_derivative1, acc.approx_derivative2, 
                             acc.approx_derivative3, acc.approx_derivative4, 
                             acc.approx_derivative5]
MainFunctions = [fcs.testing_function1, fcs.testing_function2,
                 fcs.testing_function3, fcs.testing_function4,
                 fcs.testing_function5]
MainDerivativeFunctions = [fcs.derivative_function1, fcs.derivative_function2,
                           fcs.derivative_function3, fcs.derivative_function4,
                           fcs.derivative_function5]

def get_delta_from_pow(pow):
    return 2 ** (1 - pow)

def get_dependency(approx_derivative, derivative, testing_function):
    x = []
    y = []
    for i in range(MaxPowerOfDelta):
        x.append(math.log(get_delta_from_pow(i + 1)))
        y.append(math.log(abs(approx_derivative(ConstantX, get_delta_from_pow(i + 1), testing_function) - 
                    derivative(ConstantX))))
    return x, y

def get_depend_n_plt_it(approx_derivative, derivative, testing_function, iteration):
    x, y = get_dependency(approx_derivative, derivative, testing_function)
    ax.plot(x, y, 'D-', label = "approx function " + str(iteration + 1))

print(get_delta_from_pow(MaxPowerOfDelta))
for i in range(len(MainFunctions)):
    sns.set(style='darkgrid')
    fig, ax = plt.subplots(figsize=(15, 15))
    main_func = MainFunctions[i]
    for j in range(len(ApproxDerivativeFunctions)):
        get_depend_n_plt_it(ApproxDerivativeFunctions[j], MainDerivativeFunctions[i], MainFunctions[i], j)
    #plt.title(label = "main function " + str(i + 1), fontsize=40, loc='right')
    ax.legend()
    ax.minorticks_on()
    ax.grid(which='minor', linestyle='-')
    ax.set_title("look at function " + str(i + 1), fontsize=20)
    ax.set_xlabel('log(h)', fontsize=20)
    ax.set_ylabel('log(abs)', fontsize=20)
    plt.show()
