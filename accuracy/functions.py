import math
from decimal import Decimal

def testing_function1(arg):
    return math.sin(arg ** 2)

def derivative_function1(arg):
    return 2 * arg * math.cos(arg ** 2)

def testing_function2(arg):
    return math.cos(math.sin(arg))

def derivative_function2(arg):
    return -math.cos(arg) * math.sin(math.sin(arg))

def testing_function3(arg):
    return math.exp(math.sin(math.cos(arg)))

def derivative_function3(arg):
    return testing_function3(arg) * (- math.sin(arg) * math.cos(math.cos(arg)))

def testing_function4(arg):
    return math.log(arg + 3)

def derivative_function4(arg):
    return 1 / (arg + 3)

def testing_function5(arg):
    return (arg + 3) ** 0.5

def derivative_function5(arg):
    return 1 / (2 * testing_function5(arg))
