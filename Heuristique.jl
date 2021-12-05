using JuMP
using CPLEX

function heuristique(pathFileData)
    """
    Paramètres : 
    pathFileData le chemin du fichier de donnée
    """

    nbIteMax = 100

    # Récupération des données
    data = Read_file(pathFileData)
    n = data["n"] # n le nombre de revendeurs (en comptant le fournisseur)
    l = data["l"] # l horizon de planification
    f = data["f"] # f coût fixe par période de production
    u = data["u"] # u coût unitaire
    d = data["d"] # d tableau de dimension n*l, d[i,t] exprime la demande du revendeur i au temps t
    h = data["h"] # h tableau de dimension n, h[i] exprime le coût de stockage unitaire du revendeur i
    L = data["L"] # L tableau de dimension n, L[i] exprime la capacité de stockage maximale chez le revendeur i
    M = data["C"] # M constante big M qui se doit d'être supérieure à toute valeur raisonnable que peut prendre la quantité produite sur une période
    c = matrix_cout(data) # c tableau de dimension n*n, c[i][j] exprime le coût de transport du nœud i au nœud j
    
    # PB +1 aux indices de c + ceux de LSP (h, L et L0)

    # Initialisation des paramètres SC : 
    SC = Array{Float64}(undef,n,l)
    for i=1:n
        for t=1:l
            SC[i,t] = c[1,i] + c[i,1]
        end
    end

    nbIte = 0
    while nbIte <= nbIteMax
        # Résolution du problème LSP avec prise en compte des coûts de visite
        p, y, I, q = PL_LSP("PRP_instances/A_014_ABS1_15_2.prp", SC, false)
        println("p=", p)
        println("y=", y)
        println("I=", I)
        println("q=", q)

        # Résolution, pour chaque pas de temps, du VRP
        for t=1:l
            if y[t] != 0 # Si une production est lancée à la période t, sinon ne rien faire à ce temps
            # Résoudre VRP avec les paramètres : p[t] quantité produite à la période t, 
            # y[t] vaut 1 si une production est lancée à la période t, et 0 sinon,
            # I[:,t] la quantité en stock à la fin de la période t pour chaque revendeur i, et
            # q[:,t] la quantité produite pour chaque revendeur i à la période t
            end
        end

        # Mise à jour des coûts SC

        if nbIte == 0
            bestSol = # Mémorisation de la meilleure solution
        elseif CoutbestSol < CoutSolcourante
            bestSol = Solcourante # Mémorisation de la meilleure solution
        end
        nbIte += 1
    end
    return bestSol
end

# Test 
heuristique("PRP_instances/A_014_ABS1_15_2.prp")
