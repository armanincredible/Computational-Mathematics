import math

def approx_derivative1(x, delta, function):
    return (function(x + delta) - function(x)) / delta

def approx_derivative2(x, delta, function):
    return (function(x) - function(x - delta)) / delta

def approx_derivative3(x, delta, function):
    return (function(x + delta) - function(x - delta)) / (2 * delta)

def approx_derivative4(x, delta, function):
    return (4 * (function(x + delta) - function(x - delta)) / ( 3 * 2 * delta) - 
            (function(x + 2 * delta) - function(x - 2 * delta)) / ( 3 * 4 * delta))

def approx_derivative5(x, delta, function):
    return (3 * (function(x + delta) - function(x - delta)) / ( 3 * 2 * delta) - 
            3 * (function(x + 2 * delta) - function(x - 2 * delta)) / ( 5 * 4 * delta) + 
            (function(x + 3 * delta) - function(x - 3 * delta)) / ( 10 * 6 * delta))

