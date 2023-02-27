import numpy as np
from snl_newton import newton
import math

## Determination des donnees du probleme
# Le nom generique permet de definir le nom du script
# des donnees : "(nomgen).m", le nom du fichier de resultat :
# "(nomgen).res", ainsi que le nom de l'image : "(nomgen).png"
# En retour de ce script les valeurs suivantes sont definies :

def importation():
    import donnees1

    global FFX
    global nbrsol, sol, distsol, itermax, epsF, Xmax, xlim, ylim, nbrlig, nbrcol, maxcol, maxmap, cmap

    # FFX : pointeur (handle) vers la fonction systeme traitee
    FFX = donnees1.function

    # nbrsol : nombre de solution du systeme
    nbrsol = donnees1.nbrsol

    # sol(2,nbrsol) : matrice des solutions du systeme
    sol = donnees1.sol

    # distsol : distance de discrimination entre solution
    distsol = donnees1.distsol

    # itermax : nombre max. d'iteration
    itermax = donnees1.itermax

    # epsf : valeur pour le test de conv. |F(X)| < epsF
    epsF = donnees1.epsF

    # Xmax : valeur pour le test d'arret |X| > Xmax
    Xmax = donnees1.Xmax

    # xlim : intervalle des abscisses balayees
    xlim = donnees1.xlim

    # ylim : intervalle des ordonnees balayees
    ylim = donnees1.ylim

    # nbrlig : nombre de ligne (pix) de l'image
    nbrlig = donnees1.nbrlig

    # nbrcol : nombre de colonne (pix) de l'image
    nbrcol = donnees1.nbrcol

    # maxcol : vecteur des nombre max de couleurs par solution (+1 pour la div)
    maxcol = donnees1.maxcol

    # maxmap : nbre total de couleur =Sum(maxcol)
    maxmap = donnees1.maxmap

    # cmap : palette RVB des maxmap couleurs
    cmap = donnees1.cmap


%% Chargement de la palette cmap
colormap(cmap);


## Fonction d'indexation des couleurs :

def colorMapping(num_sol, iter, itermax, maxcol, maxmap):
    math.floor((((iter-1)/itermax)**(1/7))*maxcol[num_sol]+1)+maxmap-np.sum(maxcol[num_sol+1:])

# Le principe de icol fonctionne comme suit:
#    - 0 <= (iter-1)/itermax < 1 donc apres multiplication par la
#     maxcol(numsol), le nb de nuance associé a la solution 
#    - 1 <= (iter-1)/itermax*maxcol(num_sol)+1 < nb_nuances +1
#    - Donc en passant a la partie entière on se retrouve avec un entier
#     compris entre 1 et nb_nuance. On a donc directement de numero de la
#     nuance associé a une solution.
#    - Il suffit donc de rajouter les couleurs des autres solutions afin 
#     d'avoir la bonne couleur dans la cmap 

#     Enfin, du fait que la convergence soit souvent plus rapide, et qu'on 
#     atteigne pas 100 itérations, seules les premières nuances sont utilisées.
#     Pour palier a cela, on passe a la puissance 1/7 sur (iter-1)/itermax
#     afin de ne plus prendre de plus avoir une relation linéaire entre
#     iter et nb_nuance et cela permmets ainsi de prendre plus de nuances
#     pour un nombre d'itérations faible


%% Post processing : calcul de pconv et de imoy
Mat_conv = C <= sum(maxcol(1:nbrsol)); % Matrice de booléen (1 si convergence, 0 sinon)
nb_conv = sum(sum(Mat_conv));          % On somme tous les 1 de la matrice pour avoir le nombre de cas convergent
pconv=nb_conv/(nbrlig*nbrcol);         % Convergence moyenne du système
imoy=compteur_iter/nb_conv;            % Nb d'iterations moyen (pour les cas convergents uniquement)

%% Affichage et sauvegarde de l'image
image(x,y,C); % voir doc image pour comprendre
set(gca,'dataaspectratio',[1 1 1]); % permet d'egaliser les echelles x et y
imname = strcat(nomgen,'.png'); % creation du nom du fichier image
imwrite(C,cmap,imname); % voir doc imwrite
%% Sauvegarde des resultats
% ouverture du fichier et affectation d'un identifiant :
fid = fopen(strcat(nomgen,'.txt'),'w+');
% voir doc de fprintf pour les formats %s %i %f %g ...
fprintf(fid,'image creer : %s\n',imname);
fprintf(fid,'systeme : %s\n',func2str(FFX));
fprintf(fid,'taille image : %i x %i\n',nbrlig,nbrcol);
fprintf(fid,'zone d''etude : [%f %f] x [%f %f]\n',...
xlim,ylim);
fprintf(fid,'itermax : %i\n',itermax);
fprintf(fid,'epsF : %g\n',epsF);
fprintf(fid,'Xmax : %g\n',Xmax);
fprintf(fid,'temps CPU : %.2f s\n',tfin);
fprintf(fid,'%% de conv : %.2f %%\n',floor(1e4*pconv)/100);
fprintf(fid,'nb iter moy : %.2f\n',imoy);
fclose(fid);


if __name__ == "__main__":
    importation()

    C = np.zeros() 

    ## Definition et preallocation des pixels de la zone d'image ...
    x = np.linspace(xlim[0],xlim[1],nbrcol)
    y = np.linspace(ylim[0],ylim[1],nbrlig)

    compteur_iter=0;    # Initialisation du compteur d'itérations
    t = time.time()

    for j in x:     # On parcours les lignes
        for i in y: # On parcours les colones
            X = np.array([[x],
                          [y]])    # Initialisation des coordonnées au pixel considéré
            (Xsol, conv, iter) = newton(FFX, X, epsF, itermax, Xmax) # Recherche de la solution

            if conv:             # Dans le cas de convergence
            
                for k in range (nbrsol): # Renvoi directement le numero de la solution (il y a forcément une solution)
                    if np.linalg.norm(Xsol-sol[:,k]) < distsol:
                        num_sol = k

                compteur_iter = compteur_iter + iter   # incrementation des itérations

            else:
                num_sol = nbrsol+1       # Cas de non convergence et de divergence
                
            C[j,i] = icol(num_sol,iter);  # attribution du numero de couleur de la Cmap

    elapse = time.time() - t

