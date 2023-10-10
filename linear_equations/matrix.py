import numpy as np

def norm(vec):
  return np.max(np.abs(vec))

def GetLDU(matrix):
  u = np.triu(matrix)
  l = np.tril(matrix)

  lower = matrix - u
  diagonal = l + u - matrix
  upper = matrix - l

  return lower, diagonal, upper

def MatrixElem(i: int, j: int):
  if i == 99:
    return 1
  else:
    if i == j:
      return 10
    elif i == (j + 1) or j == (i + 1):
      return 1
    else: 
      return 0

def GetMatrix():
  return np.fromfunction(np.vectorize(MatrixElem), (100, 100), dtype=np.double)

def BElem(i: int):
  return (i + 1)

def GetB():
  return np.fromfunction(np.vectorize(BElem), (100,), dtype = np.double)
