import numpy as np

# FFX est la fonction choisie
# X les coordonnées du point
# epsF valeur pour le test de |F(X)| < epsF
# Xmax valeur pour le test d'arret |X| > Xmax

def snl_newton(FFX, X, epsF, itermax, Xmax):
    iter=1
    F, FX = FFX(X)

    # On entre dans la boucle tant qu'on a pas convergé ou qu'on ne diverge pas
    while ((np.linalg.norm(F) > epsF) and (iter < itermax) and (np.linalg.norm(X) < Xmax)):
        F, FX = FFX(X)
        deltaX = -(np.linalg.inv(FX)) @ F
        X = X + deltaX
        iter+= 1

    conv = np.linalg.norm(F) <= epsF

    return (X, conv, iter)