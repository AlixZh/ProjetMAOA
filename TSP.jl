using JuMP
using CPLEX

function PL_TSP(data, client_t, affichage)
    """
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
    clients_t l'ensemble de client à traiter (on ne traite que les clients qui ont une demande pour ce pas de temps)
    affichage qui vaut true si on veut afficher les solutions trouvées
    """
    # Récupération des données
    cout = matrix_cout(data) # le cout de transport du client i (1:n+1) au j (1:n+1) (l'indice 1 est le centre de depot) 	
	# Pour rappel, tous les indices de cout doivent être augmentés de 1 par rapport aux indices de l'énoncé (en julia, les indices commencent à 1)

    N = length(client_t) # Nombre de client à traiter à ce pas de temps

    # Création d'un modèle, ce modèle fera l'interface avec le solveur CPLEX
    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_EPINT" => 1e-15 )) # Pour que les variables binaires soient bien 0 ou 1

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
    # write_to_file(m, "model_TSP_1.lp")

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

        while !is_tsp_solved(m, x)
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
    N = size(x)[1] # Va de 1 à N
    x_val = JuMP.value.(x)

    # Trouver un cycle
    cycle_idx = Int[]
    push!(cycle_idx, 1)
    while true
        v, idx = findmax(x_val[cycle_idx[end],:])
        if idx == cycle_idx[1]
            break
        else
            push!(cycle_idx, idx)
        end
    end
    
    if length(cycle_idx) < N # On ajoute la contrainte associé au cycle si il n'est pas de taille N 
        @constraint(m, sum(x[cycle_idx[i],cycle_idx[i+1]] for i in 1:length(cycle_idx)-1) + x[cycle_idx[length(cycle_idx)],cycle_idx[1]] <= length(cycle_idx)-1)
        return false
    end
    return true
end

function TSP_to_Circuit(data, client_t, x) 
    """
	Permet de représenter la solution renvoyé par le PL sous forme de liste de liste (représentant les circuits, dans l'ordre de passage)
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
    client_t le tableau de l'ensemble des noeuds (entiers) à traiter
	x les valeurs des variables binaires xij qui vallent 1 si un véhicule utilise l’arc (i, j) 
    """
	circuit = [0]

	for i in 1:length(client_t)
		# On sait que tout chemin commence par le sommet 0 (centre de dépôt)
		if x[1, i] > 0
			# On traite le circuit qui commence par l'arc (0,i)
			depart = i # On veut chercher l'arc dont l'extrémité de départ est i

			while depart != 1 # Tant qu'on ne retourne pas en 0
				# On ajoute le sommet que l'on vient de rencontrer si ce n'est pas 0 (on ne veut avoir 0 qu'au début de la liste)
				push!(circuit, client_t[depart])

				for j in 1:length(client_t)
					if j != depart
						if x[depart, j] > 0 # On a trouvé l'arc dont l'extrémité de départ est depart, c'est (départ, j)
							depart = j # C'est maintenant le départ du nouvel arc que l'on traite
							break
						end
					end
				end
			end
		end
	end

	return circuit
end

# ----- Tests -----

# pathFileData = "PRP_instances/A_014_ABS1_15_1.prp"
# data = Read_file(pathFileData)
# client_t = [0,1,2,3,4,5,6,7,8]
# x = PL_TSP(data, client_t, false)
# println("x = ", x)

# circuits = TSP_to_Circuit(data, client_t, x)
# println("circuits = ", circuits)