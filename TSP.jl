using JuMP
using CPLEX

#include("tools.jl")

function PL_TSP(data, client_t, affichage)
    """
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
    clients_t l'ensemble de client à traiter (on ne traite que les clients qui ont une demande pour ce pas de temps)
    affichage qui vaut true si on veut afficher les solutions trouvées
    """

    # Récupération des données
    # n = data["n"] # n le nombre de revendeurs (en comptant le fournisseur)
    # l = data["l"] # l horizon de planification
    # f = data["f"] # f coût fixe par période de production
    # u = data["u"] # u coût unitaire
    # d = data["d"] # d tableau de dimension n*l, d[i,t] exprime la demande du revendeur i au temps t
    # h = data["h"] # h tableau de dimension n, h[i] exprime le coût de stockage unitaire du revendeur i
    # L = data["L"] # L tableau de dimension n, L[i] exprime la capacité de stockage maximale chez le revendeur i
    # L0 = data["L0"] # L0 tableau de dimension n, L0[i] exprime la quantité en stock à la fin de la période 0 (càd au début) pour i 
    # M = data["C"] # M constante big M qui se doit d'être supérieure à toute valeur raisonnable que peut prendre la quantité produite sur une période

    cout = matrix_cout(data) # le cout de transport du client i (1:n+1) au j (1:n+1) (l'indice 1 est le centre de depot) 	
	# PB cout indice +1

    N = length(client_t) # Nombre de client à traiter à ce pas de temps

    # Création d'un modèle, ce modèle fera l'interface avec le solveur CPLEX
    m = Model(CPLEX.Optimizer)

    # Création des variables
    @variable(m, x[1:N, 1:N], Bin)

    # Fonction objectif
    @objective(m, Min, sum(x[i,j]*cout[client_t[i]+1, client_t[j]+1] for i in 1:N, j in 1:N))

    # Ajout des contraintes dans le modèle
    for i in 1:N
        @constraint(m, x[i,i] == 0)
        @constraint(m, sum(x[i, :]) == 1)
        @constraint(m, sum(x[:, i]) == 1)
        for j in 1:N
            @constraint(m, x[i,j] + x[j,i] <= 1)
        end
    end
    
    # Pour le debug :
    # Ecrit sur disque le PL au format lp
    # write_to_file(m, "model_TSP.lp")

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
    
        while !is_tsp_solved(m)
            optimize!(m)
        end

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
            
			println("\t x = ", value.(x))

            println("Temps de résolution :", solve_time(m))
        else
            println("Problème lors de la résolution")
        end
    else # Sans l'affichage, il faut qd même optimiser
        optimize!(m)
        while !is_tsp_solved(m, value.(x))
            optimize!(m)
        end
    end
    return value.(x) # On récupère les valeurs des variables de décision
end

function is_tsp_solved(m, x)
    """
    Paramètres : 
    m le modèle d TSP faisant l'interface avec le solveur CPLEX
    x l'ensemble de variables binaires indiquant si on prend l'arête ou non
    """
    N = size(x)[1]
    x_val = JuMP.value.(x)
    
    # find cycle
    cycle_idx = Int[]
    push!(cycle_idx, 1)
    while true
        v, idx = findmax(x_val[cycle_idx[end],1:N])
        if idx == cycle_idx[1]
            break
        else
            push!(cycle_idx,idx)
        end
    end
    println("cycle_idx: ", cycle_idx)
    println("Length: ", length(cycle_idx))
    if length(cycle_idx) < N
        @constraint(m, sum(x[cycle_idx,cycle_idx]) <= length(cycle_idx)-1)
        return false
    end
    return true
end

# Tests
pathFileData = "PRP_instances/A_014_ABS1_15_1.prp"
data = Read_file(pathFileData)
x = PL_TSP(data, [0,1,2], false)
println("x=", x)