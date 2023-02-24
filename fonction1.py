# fonction renvoyant la fonction vectorielle F et la
# matrice Jacobienne FX du systeme qui s'ecrit en complexe :
# z^3-1 = 0

import numpy as np

def function (X):
    x = X(1)
    y = X(2)
    F = np.zeros((1,2))
    FX = np.zeros((2,2))

    x2 = x**2
    y2 = y**2

    ## Le systeme
    F[0,0] = x * (x2 - 3*y2) - 1
    F[1,0] = y * (3*x2 - y2)

    # La matrice jacobienne
    FX[0,0] = 3 * (x2 - y2)
    FX[0,1] = -6 * x * y
    FX[1,0] = -FX[0,1]
    FX[1,1] = FX[0,0]

    return (F,FX)