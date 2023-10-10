import numpy as np
import pprint
import scipy
import scipy.linalg 

def Gauss(A,b):
    #A = our intial matrix
    #b = solution vector
    n = len(A)
    M = A

    M = np.hstack((M,np.array([b]).T))

    for i in range(n):

        leading = i + np.argmax(np.abs(A[:,i][i:]))
        M[[i, leading]] = M[[leading, i]] 

        M[i] /= M[i][i]
        row = M[i]

        for r in M[i + 1:]:
            r -= r[i] * row

    for i in range(n - 1, 0, -1):
        row = M[i]
        for r in reversed(M[:i]):
            r -= r[i] * row

    return M[:,-1]

def LUmethod(A, f):
    P, L, U = scipy.linalg.lu(A)
    y = Gauss(L, f)
    x = Gauss(U, y)
    return x
