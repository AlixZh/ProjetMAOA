using JuMP
using CPLEX

#include("tools.jl")

# const OPTIMAL = JuMP.MathOptInterface.OPTIMAL
# const INFEASIBLE = JuMP.MathOptInterface.INFEASIBLE
# const UNBOUNDED = JuMP.MathOptInterface.DUAL_INFEASIBLE;

function PL_VRP(data, demande, t, affichage)
    """
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	demande un tableau [1:n, 1:l] de taille n*l, demande[i,t] est la quantité produite pour le revendeur i à la période t, soit q[i,t] dans la notation LSP
	t l'instant de temps à considérer
	affichage qui vaut true si on veut afficher les solutions trouvées
    """

	# Récupération des données
    k = data["k"] # k le nombre de véhicules
    Q = data["Q"] # Q la capacité maximale de chaque véhicule
	n = data["n"] # n le nombre de clients

	cout = matrix_cout(data) # le cout de transport du client i (0:n) au j (0:n) (l'indice 1 est le centre de depot) 	

	# Création d'un modèle. Ce modèle fera l'interface avec le solveur CPLEX
	m = Model(CPLEX.Optimizer)

	# Création des variables
	@variable(m, x[0:n, 0:n], Bin) # x[i,j] associé à l'arête (i,j)
	# for i in 0:n # PB
    #     delete(m, x[i, i]) # On peut supprimer les variables x[i,i] comme il n'y a pas d'arête d'un noeud vers lui-même
    # end
	@variable(m, 0<= w[1:n] <=Q) # indice 1 correspond au client 1
	
	# Fonction objectif
	@objective(m, Min, sum(cout[i+1,j+1]*x[i,j] for i in 0:n for j in 0:n)) # PB soit sum(x .*cout) # C'est-à-dire 
	
	# Ajout des contraintes dans le modèle
	@constraint(m, sum(x[0,:])<=k) # Contrainte (6)
	@constraint(m, sum(x[:,0])<=k) # Contrainte (7)

	for i in 1:n
		@constraint(m, sum(x[i,:])==1) # Contrainte (8)
		@constraint(m, sum(x[:,i])==1) # Contrainte (9)

		for j in 1:n
			@constraint(m, (demande[i,t]-(Q + demande[i,t])*(1 - x[i,j]))<= (w[i]-w[j])) # Contrainte (10) MTZ
		end
	end

	# Pour le debug :
    # Ecrit sur disque le PL au format lp
    # write_to_file(m, "model_VRP.lp")

	if affichage
		# Affichages
		#for e in edges(G)  
		#if(src(e)!=1 && dst(e)!=1)
		#Affichage du modèle
		println("Affichage du modèle avant résolution:")
		print(m)
		println()

		#Résolution du problème d'optimisation linéaire m par le solveur CPLEX
		println("Résolution par le solveur linéaire choisi")
		optimize!(m)
		println()

		# Affiche tous les détails d'une solution à l'écran
		#println("Affichage de tous les détails de la solution avec la commande solution_summary")
		#println(solution_summary(m, verbose=true))
		#println()

		# Mais on peut vouloir récupérer juste une information précise
		println("Récupération et affichage \"à la main\" d'informations précises")
		status = termination_status(m)

		if status == INFEASIBLE
			println("Le problème n'est pas réalisable")
		# elseif status == UNBOUNDED
		# 	println("Le problème est non borné")
		elseif status == OPTIMAL
			println("Valeur optimale = ", objective_value(m))
			println("Solution optimale :")
			println("\t x = ", value.(x))
			println("\t x = ", value.(w))
			println("Temps de résolution :", solve_time(m))
		else
			println("Problème lors de la résolution")
		end
	else # Sans l'affichage, il faut qd même optimiser
        optimize!(m)
    end
	
	return value.(x), value.(w) # On récupère les valeurs des variables de décision
end

function VRP_to_Circuit(data, x) 
    """
	Permet de représenter la solution renvoyé par le PL sous forme de liste de liste (représentant les circuits, dans l'ordre de passage)
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	x les valeurs des variables binaires xij qui vallent 1 si un véhicule utilise l’arc (i, j) 
	filename le nom sous lequel on veut enregistrer le graphe (0 si on ne veut pas l'enregistrer) 
    """

	# Récupération des données
	n = data["n"] # n le nombre de clients	

	circuits = Vector{Int64}[]

	for i in 1:n
		# On sait que tout chemin commence par le sommet 0 (centre de dépôt)
		if x[0,i] > 0
			circuit = [0] # On traite le circuit qui commence par l'arc (0,i)
			depart = i # On veut chercher l'arc dont l'extrémité de départ est i

			while depart != 0 # Tant qu'on ne retourne pas en 0
				# On ajoute le sommet que l'on vient de rencontrer si ce n'est pas 0 (on ne veut avoir 0 qu'au début de la liste)
				push!(circuit, depart)

				for j in 0:n
					if j != depart
						if x[depart, j] > 0 # On a trouvé l'arc dont l'extrémité de départ est depart, c'est (départ, j)
							depart = j # C'est maintenant le départ du nouvel arc que l'on traite
							break
						end

					end
				end

			end

			push!(circuits, circuit)

		end
	end

	# if filename != 0
	# 	#PB  à tester
	# 	# On enregistre le graphe dans le répertoire courant 
	# 	WritePngGraph_Boites(data, q, t, circuits, filename)
	# end

	return circuits

end

# # Tests
# pathFileData = "PRP_instances/A_014_ABS1_15_1.prp"
# data = Read_file(pathFileData)

# # Récupérer le q dans LSP
# p, y, I, q = PL_LSP(data, 0, false)

# x, w = PL_VRP(data, q, 2, false) # Le q vient de LSP.jl
# println("x=", x)
# println("w=", w)

# circuits = VRP_to_Circuit(data, q, x)
# println("Circuits =", circuits)

# PB Voir VRP_Heuristique pour évaluer le coût des circuits

# #PB tests
# boites = [[0,1,2],[4,5,6],[3,7,8]]
# temp_boites = boites
# for ind_boite in 1:length(boites)
# 	boites[ind_boite] = []
# end
# println("boites = ", boites)
# println("temp_boites = ", temp_boites)



# bo = [1,2,0,1,2,7,0,1,3]
# # filter!(e->e!=0, bo)
# insert!(bo, 1, 9)
# println(bo)
println(0%3)
println(1%3)
println(2%3)
println(3%3)