import numpy as np
import matplotlib.pyplot as plt

#XI.9.5

def f(x):
    return np.cos(2 * np.pi * x)

def P2(x):
    return 10 + np.sin(2 * np.pi * x)

def progon(A, r, N, h):
    h2 = h**2
    for n in range(1, N - 1):
        A[n][n - 1] = 1.0/h2
        A[n][n]     = -2.0/h2 - P2(n * h)
        A[n][n + 1] = 1.0/h2

        r[n] = f(n * h)
   
    A[0][0]     = -2.0/h2 - P2(0)
    A[0][1]     = 1.0/h2
    A[0][N - 1] = 1.0/h2
    r[0]        = f(0)

    A[N - 1][0]     = 1.0/h2
    A[N - 1][N - 2] = 1.0/h2
    A[N - 1][N - 1] = -2.0/h2 - P2(1.0 - h)
    r[N - 1] = f(1.0 - h)


def main():
    h = 0.005
    N = int(1.0/h)
    assert (N * h == 1)
    A = np.zeros((N, N))
    r = np.zeros(N)

    progon(A, r, N, h)

    ys = np.linalg.solve(A, r)
    xs = np.arange(0, 1.0, h)

    plt.grid()
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Final solution without periods")
    plt.plot(xs, ys)
    plt.show()

    repetitions = 5

    plt.grid()
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title(f"Final solution with {repetitions} periods")
    plt.plot(np.arange(0, 1.0 * repetitions, h), np.tile(ys, repetitions))
    plt.show()


if __name__ == "__main__":
  main()
