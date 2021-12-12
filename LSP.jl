using JuMP
using CPLEX

#include("tools.jl")

function PL_LSP(data, SC, affichage)
    """
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
    SC sert pour l'heuristique du PDI et représente le coût de visite du client i à la période t (vaut 0 si on n'est pas dans le cas heuristique du PDI)
    affichage qui vaut true si on veut afficher les solutions trouvées
    """

    # Récupération des données
    n = data["n"] # n le nombre de revendeurs (en comptant le fournisseur)
    l = data["l"] # l horizon de planification
    f = data["f"] # f coût fixe par période de production
    u = data["u"] # u coût unitaire
    k = data["k"] # k le nombre de véhicules
    Q = data["Q"] # Q la capacité maximale de chaque véhicule
    d = data["d"] # d tableau de dimension n*l, d[i,t] exprime la demande du revendeur i au temps t
    h = data["h"] # h tableau de dimension n, h[i] exprime le coût de stockage unitaire du revendeur i
    L = data["L"] # L tableau de dimension n, L[i] exprime la capacité de stockage maximale chez le revendeur i
    L0 = data["L0"] # L0 tableau de dimension n, L0[i] exprime la quantité en stock à la fin de la période 0 (càd au début) pour i 
    M = data["C"] # M constante big M qui se doit d'être supérieure à toute valeur raisonnable que peut prendre la quantité produite sur une période
    # Pour rappel, tous les indices de h, L et L0 doivent être augmentés de 1 par rapport aux indices de l'énoncé (en julia, les indices commencent à 1)

    # Création d'un modèle, ce modèle fera l'interface avec le solveur CPLEX
    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_EPINT" => 1e-15 )) # Pour que les variables binaires soient bien 0 ou 1

    # Création des variables
    @variable(m, 0<=p[1:l]) # p tableau de dimension l, p[t] exprime la quantité produite à la période t
    @variable(m, y[1:l], Bin) # y tableau de dimension l, y[t] vaut 1 si une production est lancée à la période t, et 0 sinon
    @variable(m, 0<=I[0:n, 0:l]) # I tableau de dimension (n+1)*(l+1), I[i,j] exprime la quantité en stock à la fin de la période j pour le revendeur i
    @variable(m, 0<=q[1:n, 1:l]) # q tableau de dimension n*l, q[i,j] exprime la quantité produite pour le revendeur i à la période j
    if SC != 0
        @variable(m, z[1:n, 1:l], Bin) 
    end

    # Fonction objectif
    if SC == 0
        # LSP normal (sans heuristique PDI)
        @objective(m, Min, sum(u*p[t] + f*y[t] + sum(h[i+1]*I[i,t] for i = 1:n) for t = 1:l ) )
    else
        # LSP dans le cadre de l'heuristique PDI
        @objective(m, Min, sum(u*p[t] + f*y[t] + sum(h[i+1]*I[i,t]+SC[i,t]*z[i,t] for i = 1:n) for t = 1:l ) )
        
        # Ajout des contraintes couplantes entre les variables utiles à l'heuristique PDI et les variables de production dans le modèle
        for i in 1:n
            for t in 1:l
                @constraint(m, q[i,t]<=z[i,t]*M) # Contrainte big M qui exprime que z[i,t] vaut 1 si le client i est visité à la période t 
                                                 # (càd que l'on a produit pour i à la période t) comme c'est en minimisation, 
                                                 # z[i,t] vaudra 0 si q[i,t] vaut 0, on a donc bien ce que l'on voulait
                
                # Pour être sur que si on ne produit pas, z est nul, on peut ajouter la contrainte suivante :
                # @constraint(m, z[i,t] <= q[i,t])
            end
        end
    end
   
   # Ajout des contraintes dans le modèle
   for i in 0:n
       @constraint(m, I[i, 0] == L0[i+1]) # Initialisation des I au temps 0
   end

   for t in 1:l
       @constraint(m, I[0,t-1] + p[t] == sum(q[i,t] for i = 1:n) + I[0,t] ) # Contrainte (1)
       @constraint(m, p[t] <= M*y[t]) # Contrainte (3)
       @constraint(m, I[0,t-1] <= L[1]) # Contrainte (4) 
       @constraint(m, sum(q[i,t] for i in 1:n) <= Q*k) # Ajout de cette contrainte pour que le VRP soit réalisable (la quantité à distribuer ne doit pas excéder la capacité maximale totale pour chaque pas de temps)

       for i in 1:n
           @constraint(m, I[i,t-1] + q[i,t] == d[i,t] + I[i,t] ) # Contrainte (2)
           @constraint(m, I[i,t-1] + q[i,t] <= L[i+1] ) # Contrainte (5)
        end
    # On n'en n'a pas besoin ici mais il peut être important de nommer les variables pour récupérer la solution duale
    # Pour l'afficher, on fera println("\t c1 = ", -dual(c1)) # la fonction dual renvoie le coût réduit qui est l'opposé de la variable duale en maximisation
    end

    # Pour le debug :
    # Ecrit sur disque le PL au format lp
    # write_to_file(m, "model_LSP.lp")

    if affichage
    # Affichages
    # Affichage du modèle
        # println("Affichage du modèle avant résolution:")
        # print(m)
        # println()

        # Résolution du problème d'optimisation linéaire m par le solveur CPLEX
        println("Résolution par le solveur linéaire choisi")
        optimize!(m)
        println()
    
        # Affiche tous les détails d'une solution à l'écran
        println("Affichage de tous les détails de la solution avec la commande solution_summary")
        println(solution_summary(m, verbose=true))
        println()

        # Mais on peut vouloir récupérer juste une information précise
        println("Récupération et affichage \"à la main\" d'informations précises")
        status = termination_status(m)

        if status == INFEASIBLE
            println("Le problème n'est pas réalisable")
        #elseif status == UNBOUNDED
        #    println("Le problème est non borné")
        elseif status == OPTIMAL # ou JuMP.MathOptInterface.OPTIMAL
            println("Valeur optimale = ", objective_value(m))
            println("Solution optimale :")
            
			println("\t p = ", value.(p))
			println("\t y = ", value.(y))
			println("\t I = ", value.(I))
			println("\t q = ", value.(q))
            
            # for i= 1:l
            #     print("\t p[",i,"] = ", value(p[i]))
            # end
            # println("")
            # for i= 1:l
            #     print("\t y[",i,"] = ", value(y[i]))
            # end
            # println("")
            # for i= 1:n
            #     for j= 1:l
            #         print("\t I[",i,",",j,"] = ", value(I[i,j]))    
            #     end
            #     println("")
            # end
            # for i= 1:n
            #     for j= 1:l
            #         print("\t q[",i,",",j,"] = ", value(q[i,j]))
            #     end
            #     println("")
            # end
            println("Temps de résolution :", solve_time(m))
        else
            println("Problème lors de la résolution")
        end
    else # Sans l'affichage, il faut qd même optimiser
        optimize!(m)
    end
    return value.(p), value.(y), value.(I), value.(q) # On récupère les valeurs des variables de décision
end

# ----- Tests -----

pathFileData = "PRP_instances/A_014_ABS1_15_1.prp"
data = Read_file(pathFileData)
p, y, I, q = PL_LSP(data, 0, false)
println("p=", p)
println("y=", y)
println("I=", I)
println("q=", q)

# PB régler les import en haut des fichiers
# PB régler les commentaires dans les tests (tout mettre en commentaire sauf PDI)