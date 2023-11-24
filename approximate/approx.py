import matplotlib.pyplot as plt
import numpy as np
import bisect

class InterpolationMethod:
    def __init__(self, x_data, y_data, name):
        self.x_data = x_data
        self.y_data = y_data
        self.name = name

    def get_name(self):
        return self.name

    def get_value(self, x):
        raise NotImplementedError("Subclasses must implement get_value method")

class NewtonApprox(InterpolationMethod):
    def __init__(self, x_data, y_data):
        super().__init__(x_data, y_data, "Newton")
        self.b_coef = np.zeros((len(x_data), len(x_data)))

        for i in range(1, len(x_data)):
            for j in range(len(x_data) - i):
                if i == 1:
                    self.b_coef[i][j] = (y_data[j + 1] - y_data[j]) / (x_data[j + 1] - x_data[j])
                else:
                    self.b_coef[i][j] = (self.b_coef[i - 1][j + 1] - self.b_coef[i - 1][j]) / (x_data[j + i] - x_data[j])

    def get_value(self, x):
        result = self.y_data[0]

        for i in range(1, len(self.x_data)):
            res = self.b_coef[i][0]
            for j in range(1, i + 1):
                res *= (x - self.x_data[j - 1])
            result += res

        return result

class SplineApprox(InterpolationMethod):
    def __init__(self, x_data, y_data):
        super().__init__(x_data, y_data, "Spline")
        self.M_coefs = np.zeros(len(x_data))

        eq_matrix = np.zeros((len(x_data) - 2, len(x_data) - 2))
        for i in range(len(x_data) - 2):
            for j in (-1, 0, 1):
                if (i + j) >= (len(x_data) - 2) or (i + j) < 0:
                    continue
                eq_matrix[i][i + j] = ((x_data[i + 1] - x_data[i]) * (j != 1) + (x_data[i + 2] - x_data[i + 1]) * (j != -1)) / (3 * (1 + (j != 0)))

        rhs = [((y_data[i + 2] - y_data[i + 1]) / (x_data[i + 2] - x_data[i + 1]) - (y_data[i + 1] - y_data[i]) / (x_data[i + 1] - x_data[i])) for i in range(len(x_data) - 2)]

        coefs = np.linalg.solve(eq_matrix, rhs)
        for i in range(len(x_data) - 2):
            self.M_coefs[i + 1] = coefs[i]

    def get_value(self, x):
        pos = min(bisect.bisect_right(self.x_data, x), len(self.x_data) - 1)

        M0 = self.M_coefs[pos - 1]
        M1 = self.M_coefs[pos]
        x0 = self.x_data[pos - 1]
        x1 = self.x_data[pos]
        h = x1 - x0
        f1 = self.y_data[pos]
        f0 = self.y_data[pos - 1]

        return (x - x0) ** 3 / (6 * h) * M1 + (x1 - x) ** 3 / (6 * h) * M0 + (x - x0) / h * (f1 - h**2 * M1 / 6) + (x1 - x) / h * (f0 - h**2 * M0 / 6)

class LeastSqApprox(InterpolationMethod):
    def __init__(self, x_data, y_data):
        super().__init__(x_data, y_data, "LeastSquares")
        phi = [lambda x: x ** 2, lambda x: x, lambda x: 1]

        matrix = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                for x_value in self.x_data:
                    matrix[i][j] += phi[i](x_value) * phi[j](x_value)
        rhs = [0.0 for _ in range(3)]
        for i in range(3):
            for j in range(len(self.x_data)):
                rhs[i] += self.y_data[j] * (self.x_data[j] ** (2 - i))

        self.coefs = np.linalg.solve(matrix, rhs)

    def get_value(self, x):
        return self.coefs[0] * x**2 + self.coefs[1] * x + self.coefs[2]

def main():
    data = {
        1910: 92228496,
        1920: 106021537,
        1930: 123202624,
        1940: 132164569,
        1950: 151325798,
        1960: 179323175,
        1970: 203211926,
        1980: 226545805,
        1990: 248709873,
        2000: 281421906
    }

    y_data = list(data.values())
    x_data = list(data.keys())

    methods = [NewtonApprox(x_data, y_data), SplineApprox(x_data, y_data), LeastSqApprox(x_data, y_data)]
    x_data_plot = list(range(1910, 2020, 10))

    for method in methods:
        new_y_data = [method.get_value(x_value) for x_value in x_data_plot]

        plt.title(method.get_name())
        plt.ylabel("Y")
        plt.xlabel("X")
        plt.grid()

        plt.plot(x_data_plot, new_y_data)
        plt.savefig("img/" + method.get_name())
        plt.show()

if __name__ == '__main__':
    main()
