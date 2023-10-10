import numpy as np
import matrix as mt

def Jacobi(A, b, N=40, x=None):
    #A = our initial matrix
    #b = solution vector
    #N = number of iterations
    # Create an initial guess of vector 0                                                                                                                                                           
    if x is None:
        x = np.zeros(len(A[0]))

    data = np.zeros(N)

    # Create a vector of the np.diagonal elements of A                                                                                                                                                
    # and subtract them from A                                                                                                                                                                     
    D = np.diag(A)
    T = A - np.diagflat(D)
    
    # Iterate for N times                                                                                                                                                                          
    for i in range(N):
        x = (b - np.dot(T,x)) / D   
        diff = b - np.matmul(A, x)
        data[i] = mt.norm(diff)  

    print("The number of iteration is: %d" %(i+1))
    return x, data


def GaussSeidel(A, b, N=25, x=None):
    
    if x is None:
        x = np.zeros(len(A[0]))

    data = np.zeros(N)

    L = np.tril(A)
    U = A - L
    for i in range(N):
        x = np.dot(np.linalg.inv(L), b - np.dot(U, x))
        diff = b - np.matmul(A, x)
        data[i] = mt.norm(diff)

    return x, data

def UpperRelaxation(A, b, w, N=20, x=None):
    
    if x is None:
        x = np.zeros(len(A[0]))

    data = np.zeros(N)

    u = np.triu(A)
    l = np.tril(A)
    L = A - u
    D = l + u - A
    U = A - l

    B = - np.matmul(np.linalg.inv(D + w * L), (w - 1) * D + w * U)
    F =   np.matmul(np.linalg.inv(D + w * L), b) * w

    for i in range(N):
        x = np.matmul(B, x) + F
        diff = b - np.matmul(A, x)
        data[i] = mt.norm(diff)

    return x, data