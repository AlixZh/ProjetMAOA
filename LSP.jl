using JuMP
using CPLEX

function PL_LSP(pathFileData, SC, affichage)
    """
    Paramètres : 
    pathFileData le chemin du fichier de donnée
    SC sert pour l'heuristique et représente le coût de visite du client i à la période t (vaut 0 si on n'est pas dans le cas heuristique)
    affichage qui vaut true si on veut afficher les solutions trouvées
    """

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

#   L0 = data["L0"] pour cste en + ?

    #PB +1 aux indices de h, L et L0 ?
#    println("h=",h)
#    println("h(1)=",h[1])
#    println("l=",l)
#    println("n=",n)
#    println("L=",L)
#    println("L[1]=",L[1])
#    println("data[L0]=",data["L0"])
#    println("h(n)=",h[n])
#    println("u=",u)
#    println("f=",f)

    # Création d'un modèle, ce modèle fera l'interface avec le solveur CPLEX
    m = Model(CPLEX.Optimizer)

    # Création des variables
    @variable(m, 0<=p[1:l]) # p tableau de dimension l, p[t] exprime la quantité produite à la période t
    @variable(m, y[1:l], Bin) # y tableau de dimension l, y[t] vaut 1 si une production est lancée à la période t, et 0 sinon
    @variable(m, 0<=I[0:n, 0:l]) # I tableau de dimension (n+1)*(l+1), I[i,j] exprime la quantité en stock à la fin de la période j pour le revendeur i
    @variable(m, 0<=q[1:n, 1:l]) # q tableau de dimension n*l, q[i,j] exprime la quantité produite pour le revendeur i à la période j
    if SC != 0
        @variable(m, z[1:n, 1:l], Bin) # PB ajouter des contraintes couplantes entre ces variables et les variables de production dans le modele
    end

    # Fonction objectif
    if SC == 0
        println("PB SANS heuristique")
        @objective(m, Min, sum(u*p[t] + f*y[t] + sum(h[i+1]*I[i,t] for i = 1:n) for t = 1:l ) )
    else
        println("PB AVEC heuristique")
        @objective(m, Min, sum(u*p[t] + f*y[t] + sum(h[i+1]*I[i,t]+SC[i,t]*z[i,t] for i = 1:n) for t = 1:l ) )
    end
   
   # Ajout des contraintes dans le modèle

   for t in 1:l
       @constraint(m, I[0,t-1] + p[t] == sum(q[i,t] for i = 1:n) + I[0,t] ) # Contrainte (1)
       @constraint(m, p[t] <= M*y[t]) # Contrainte (3)
       @constraint(m, I[0,t-1] <= L[1]) # Contrainte (4) 

       for i in 1:n
           @constraint(m, I[i,t-1] + q[i,t] == d[i,t] + I[i,t] ) # Contrainte (2)
           @constraint(m, I[i,t-1] + q[i,t] <= L[i+1] ) # Contrainte (5)
        end
    # On n'en n'a pas besoin ici mais il peut être important de nommer les variables pour récupérer la solution duale
    # Pour l'afficher, on fera println("\t c1 = ", -dual(c1)) # la fonction dual renvoie le coût réduit qui est l'opposé de la variable duale en maximisation
    end

    if affichage
    # Affichages
    # Affichage du modèle
        # println("Affichage du modèle avant résolution:")
        # print(m)
        # println()

        # Résolution du problème d'optimisation linéaire m par le solveur GLPK
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
            for i= 1:l
                print("\t p[",i,"] = ", value(p[i]))
            end
            println("")
            for i= 1:l
                print("\t y[",i,"] = ", value(y[i]))
            end
            println("")
            for i= 1:n
                for j= 1:l
                    print("\t I[",i,",",j,"] = ", value(I[i,j]))    
                end
                println("")
            end
            for i= 1:n
                for j= 1:l
                    print("\t q[",i,",",j,"] = ", value(q[i,j]))
                end
                println("")
            end
            println("Temps de résolution :", solve_time(m))
        else
            println("Problème lors de la résolution")
        end
    else # Sans l'affichage, il faut qd même optimiser
        optimize!(m)
    end

    return value.(p), value.(y), value.(I), value.(q) # On récupère les valeurs des variables de décision
end

# Test
p, y, I, q = PL_LSP("PRP_instances/A_014_ABS1_15_2.prp", 0, false)
println("p=", p)
println("y=", y)
println("I=", I)
println("q=", q)

#PB y vaut 1 toujours mm qd p vaut 0
#PB p vaut tjr 0 ??