import numpy as np
import math

# FFX : pointeur (handle) vers la fonction système traitée
# fonction renvoyant la fonction vectorielle F et la
# matrice Jacobienne FX du systeme qui s'ecrit en complexe :
# z^3-1 = 0
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

# nbrsol : nombre de solution du système
nbrsol = 3

# sol(2,nbrsol) : matrice des solutions du système
sol=np.array[[1,-0.5,-0.5],\
             [0,math.sqrt(3)/2,-math.sqrt(3)/2]]

# distsol : distance de discrimination entre solutions
distsol=np.linalg.norm(sol[:,1]-sol[:,2])/2

# itermax : nombre max. d'itération
itermax=100

# epsf : valeur pour le test de conv. |F(X)| < epsf
epsF=0.00001

# Xmax : valeur pour le test d'arrêt |X| > Xmax
Xmax=1000

# xlim : intervalle des abscisses balayées
xlim=[-2,2]

# ylim : intervalle des ordonnées balayées
ylim=[-2,2]

# nbrlig : nombre de ligne (pix) de l'image
nbrlig=1000

# nbrcol : nombre de colonne (pix) de l'image
nbrcol=1000

# maxcol : vecteur des nombre max de couleurs par solution (+1 pour la div)
maxcol = np.array([40, 40, 40, 1]).reshape(1,4) # 1 couleur par solution + 1 couleur pour la divergence

# maxmap : nombre total de couleur = sum(maxcol)
maxmap = sum(maxcol) # nombre total de couleurs

# cmap : palette RVB des maxmap couleurs
cmap =[1.0000    1.0000    1.0000
    0.9744    0.9744    0.9846
    0.9487    0.9487    0.9692
    0.9231    0.9231    0.9538
    0.8974    0.8974    0.9385
    0.8718    0.8718    0.9231
    0.8462    0.8462    0.9077
    0.8205    0.8205    0.8923
    0.7949    0.7949    0.8769
    0.7692    0.7692    0.8615
    0.7436    0.7436    0.8462
    0.7179    0.7179    0.8308
    0.6923    0.6923    0.8154
    0.6667    0.6667    0.8000
    0.6410    0.6410    0.7846
    0.6154    0.6154    0.7692
    0.5897    0.5897    0.7538
    0.5641    0.5641    0.7385
    0.5385    0.5385    0.7231
    0.5128    0.5128    0.7077
    0.4872    0.4872    0.6923
    0.4615    0.4615    0.6769
    0.4359    0.4359    0.6615
    0.4103    0.4103    0.6462
    0.3846    0.3846    0.6308
    0.3590    0.3590    0.6154
    0.3333    0.3333    0.6000
    0.3077    0.3077    0.5846
    0.2821    0.2821    0.5692
    0.2564    0.2564    0.5538
    0.2308    0.2308    0.5385
    0.2051    0.2051    0.5231
    0.1795    0.1795    0.5077
    0.1538    0.1538    0.4923
    0.1282    0.1282    0.4769
    0.1026    0.1026    0.4615
    0.0769    0.0769    0.4462
    0.0513    0.0513    0.4308
    0.0256    0.0256    0.4154
         0         0    0.4000
    1.0000    1.0000    1.0000
    0.9897    0.9744    0.9744
    0.9795    0.9487    0.9487
    0.9692    0.9231    0.9231
    0.9590    0.8974    0.8974
    0.9487    0.8718    0.8718
    0.9385    0.8462    0.8462
    0.9282    0.8205    0.8205
    0.9179    0.7949    0.7949
    0.9077    0.7692    0.7692
    0.8974    0.7436    0.7436
    0.8872    0.7179    0.7179
    0.8769    0.6923    0.6923
    0.8667    0.6667    0.6667
    0.8564    0.6410    0.6410
    0.8462    0.6154    0.6154
    0.8359    0.5897    0.5897
    0.8256    0.5641    0.5641
    0.8154    0.5385    0.5385
    0.8051    0.5128    0.5128
    0.7949    0.4872    0.4872
    0.7846    0.4615    0.4615
    0.7744    0.4359    0.4359
    0.7641    0.4103    0.4103
    0.7538    0.3846    0.3846
    0.7436    0.3590    0.3590
    0.7333    0.3333    0.3333
    0.7231    0.3077    0.3077
    0.7128    0.2821    0.2821
    0.7026    0.2564    0.2564
    0.6923    0.2308    0.2308
    0.6821    0.2051    0.2051
    0.6718    0.1795    0.1795
    0.6615    0.1538    0.1538
    0.6513    0.1282    0.1282
    0.6410    0.1026    0.1026
    0.6308    0.0769    0.0769
    0.6205    0.0513    0.0513
    0.6103    0.0256    0.0256
    0.6000         0         0
    1.0000    1.0000    1.0000
    0.9744    0.9846    0.9795
    0.9487    0.9692    0.9590
    0.9231    0.9538    0.9385
    0.8974    0.9385    0.9179
    0.8718    0.9231    0.8974
    0.8462    0.9077    0.8769
    0.8205    0.8923    0.8564
    0.7949    0.8769    0.8359
    0.7692    0.8615    0.8154
    0.7436    0.8462    0.7949
    0.7179    0.8308    0.7744
    0.6923    0.8154    0.7538
    0.6667    0.8000    0.7333
    0.6410    0.7846    0.7128
    0.6154    0.7692    0.6923
    0.5897    0.7538    0.6718
    0.5641    0.7385    0.6513
    0.5385    0.7231    0.6308
    0.5128    0.7077    0.6103
    0.4872    0.6923    0.5897
    0.4615    0.6769    0.5692
    0.4359    0.6615    0.5487
    0.4103    0.6462    0.5282
    0.3846    0.6308    0.5077
    0.3590    0.6154    0.4872
    0.3333    0.6000    0.4667
    0.3077    0.5846    0.4462
    0.2821    0.5692    0.4256
    0.2564    0.5538    0.4051
    0.2308    0.5385    0.3846
    0.2051    0.5231    0.3641
    0.1795    0.5077    0.3436
    0.1538    0.4923    0.3231
    0.1282    0.4769    0.3026
    0.1026    0.4615    0.2821
    0.0769    0.4462    0.2615
    0.0513    0.4308    0.2410
    0.0256    0.4154    0.2205
         0    0.4000    0.2000
    0.0000    0.0000    0.0000]