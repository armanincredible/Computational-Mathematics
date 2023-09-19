import math
import matplotlib.pyplot as plt

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
    plt.plot(x, y, label = "approx function " + str(iteration + 1))

print(get_delta_from_pow(MaxPowerOfDelta))
for i in range(len(MainFunctions)):
    plt.figure(figsize=(15, 15))
    main_func = MainFunctions[i]
    for j in range(len(ApproxDerivativeFunctions)):
        get_depend_n_plt_it(ApproxDerivativeFunctions[j], MainDerivativeFunctions[i], MainFunctions[i], j)
    #plt.title(label = "main function " + str(i + 1), fontsize=40, loc='right')
    plt.legend()
    plt.title("look at function " + str(i + 1))
    plt.show()
