# ProjetMAOA

## Type des données

n : nombre de clients
    int
  
g : un graphe complet non orienté, n+1 noeuds

l : nombre d'horizon

u : coût unitaire de production

f : coût fixe de production

C : capacité de production

Q : capacité maximale de chaque véhicule

k : flotte, nombre de véhicule
    int, chaque véhicule est noté par un indice de 1 à k

d : demande du client i à la période t, 1<=i<=n, 1<=t<=l
    matrice de taille n*l

L0 : capacité de stockage maximale au noeud 0 (centre de dépot) pour chacun des clients
     vecteur de taille n+1, L0[i+1] : capacité de stockage maximale au noeud 0 pour le noeud i, 0 <= i <=n
     
L : capacité de stockage maximale du noeud i
     vecteur de taille n+1, Li[i+1] : capacité de stockage maximale du noeud i, 0 <= i <= n
     
h : coût de stockage unitaire
    vecteur de taille n+1, h[i] : coût de stockage unitaire du noeud i,  0 <= i <= n
    
## Fichiers 

tools.jl : fichier contenant des fonctions pour récupérer les instances, calculer les coûts de transport pour les instances A et B
