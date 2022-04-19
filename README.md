# ProjetMAOA

Projet réalisé dans le cadre de l'UE Modèles et applications en ordonnancement et optimisation combinatoire (MAOA).
Sujet du projet : Heuristiques et PLNE compacts ou à nombre exponentiel de contraintes (Branch-and-Cut)

Ce repository est constitué de plusieurs fichiers : les fichiers de codes (BC.jl, LSP.jl, PDI1.jl, PDI2.jl, PDI_Exact.jl, PDI_Heuristique.jl, TSP.jl, Test.jl, TestPDI.jl, VRP.jl, VRP_Heuristiques.jl, tools.jl, toolsPDI.jl), MAOA_PDI_Rapport.pdf le rapport final du projet, et MAOA_Sujet.pdf l'énoncé du projet.

## Type des données

Les données des instances seront stockées dans un dictionnaire. Les clés du disctionnaire sont : n, l, u, f, C, Q, k, d, L0, L, h, coord et mc (si instance B).

n : nombre de clients
    int

l : nombre d'horizon

u : coût unitaire de production

f : coût fixe de production

C : capacité de production

Q : capacité maximale de chaque véhicule

k : flotte, nombre de véhicule
    int, chaque véhicule est noté par un indice de 1 à k

d : demande du client i à la période t, 1<=i<=n, 1<=t<=l
    matrice de taille n*l

L0 : capacité de stockage initial (t = 0) pour chacun des clients
     vecteur de taille n+1, L0[i+1] : capacité de stockage maximale au noeud 0 pour le noeud i, 0 <= i <=n
     
L : capacité de stockage maximale du noeud i
     vecteur de taille n+1, Li[i+1] : capacité de stockage maximale du noeud i, 0 <= i <= n
     
h : coût de stockage unitaire
    vecteur de taille n+1, h[i] : coût de stockage unitaire du noeud i,  0 <= i <= n
    
coord : coordonnées des positions des clients
        vecteur de tuple de taille n+1

cij : coût de transport de i à j, calculer par des fonctions à partir des coordonnées de i et j
    
## Fichiers 

tools.jl : fichier contenant des fonctions pour récupérer les instances, calculer les coûts de transport pour les instances A et B


## Lien du rapport
https://fr.overleaf.com/read/pxwfbwqkvnym 
