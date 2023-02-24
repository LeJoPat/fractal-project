clear variables
close all
%% Determination des donnees du probleme
% Le nom generique permet de definir le nom du script
% des donnees : "(nomgen).m", le nom du fichier de resultat :
% "(nomgen).res", ainsi que le nom de l'image : "(nomgen).png"
nomgen = input('Nom du cas a traiter :','s');
eval(nomgen); % execution du script "(nomgen).m"
% En retour de ce script les valeurs suivantes sont definies :
% FFX : pointeur (handle) vers la fonction systeme traitee
% nbrsol : nombre de solution du systeme
% sol(2,nbrsol) : matrice des solutions du systeme
% distsol : distance de discrimination entre solution
% itermax : nombre max. d'iteration
% epsf : valeur pour le test de conv. |F(X)| < epsF
% Xmax : valeur pour le test d'arret |X| > Xmax
% xlim : intervalle des abscisses balayees
% ylim : intervalle des ordonnees balayees
% nbrlig : nombre de ligne (pix) de l'image
% nbrcol : nombre de colonne (pix) de l'image
% maxcol : vecteur des nombre max de couleurs par solution (+1 pour la div)
% maxmap : nbre total de couleur =Sum(maxcol)
% cmap : palette RVB des maxmap couleurs

%% Chargement de la palette cmap
colormap(cmap);

%% Definition et preallocation des pixels de la zone d'image ...
x = linspace(xlim(1),xlim(2),nbrcol);
y = linspace(ylim(1),ylim(2),nbrlig);

%% Fonction d'indexation des couleurs : icol ..

icol=@(num_sol, iter) floor((((iter-1)/itermax)^(1/7))*maxcol(num_sol)+1)+maxmap-sum(maxcol(num_sol:end));

% Le principe de icol fonctionne comme suit:
%    - 0 <= (iter-1)/itermax < 1 donc apres multiplication par la
%     maxcol(numsol), le nb de nuance associé a la solution 
%    - 1 <= (iter-1)/itermax*maxcol(num_sol)+1 < nb_nuances +1
%    - Donc en passant a la partie entière on se retrouve avec un entier
%     compris entre 1 et nb_nuance. On a donc directement de numero de la
%     nuance associé a une solution.
%    - Il suffit donc de rajouter les couleurs des autres solutions afin 
%     d'avoir la bonne couleur dans la cmap 

%     Enfin, du fait que la convergence soit souvent plus rapide, et qu'on 
%     atteigne pas 100 itérations, seules les premières nuances sont utilisées.
%     Pour palier a cela, on passe a la puissance 1/7 sur (iter-1)/itermax
%     afin de ne plus prendre de plus avoir une relation linéaire entre
%     iter et nb_nuance et cela permmets ainsi de prendre plus de nuances
%     pour un nombre d'itérations faible

%% Initialisation du calcul ...
C = zeros(nbrlig, nbrcol);

%% Les boucles de calculs
compteur_iter=0;     % Initialisation du compteur d'itérations
tic % on declenche le chrono (cf. doc tic toc pour l'unité)
for j=1 :nbrlig      % On parcours les lignes
    for i=1 :nbrcol  % On parcours les colones
        X = [x(i);y(j)];     % Initialisation des coordonnées au pixel considéré
        [X,conv,iter]= SNL_NEWTON(FFX,X,epsF,itermax,Xmax);  % Recherche de la solution
        if conv == 1             % Dans le cas de convergence
            for k=1:nbrsol  
                T(k) = norm(X-sol(:,k))<distsol; % T est une liste de booléens de taille nbrsol (1 pour la solution)
            end
            num_sol = find(T);  % Renvoi directement le numero de la solution (1 apparait une seule fois dans T)
            compteur_iter = compteur_iter+iter;  % incrementation des itérations
        else 
            num_sol = nbrsol+1;       % Cas de non convergence et de divergence
        end
        C(j,i) = real(icol(num_sol,iter));  % attribution du numero de couleur de la Cmap
    end
end
            
tfin = toc; % on enregistre le temps de calcul

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


import numpy as np
C = np.zeros()


if __name__ == "__main__":

