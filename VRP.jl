using JuMP
using CPLEX

include("tools.jl")
# const OPTIMAL = JuMP.MathOptInterface.OPTIMAL
# const INFEASIBLE = JuMP.MathOptInterface.INFEASIBLE
# const UNBOUNDED = JuMP.MathOptInterface.DUAL_INFEASIBLE;

function PL_VRP(pathFileData, demande, t, affichage)
    """
    Paramètres : 
    pathFileData le chemin du fichier de donnée
	demande un tableau [1:n, 1:l] de taille n*l, demande[i,t] est la demande du client i au temps t, soit q[i,t] dans la notation LSP
	t l'instant de temps à considérer
	affichage qui vaut true si on veut afficher les solutions trouvées
    """

	# Récupération des données
    data = Read_file(pathFileData)
    k = data["k"] # k le nombre de véhicules
    Q = data["Q"] # Q la capacité maximale de chaque véhicule
	n = data["n"] # n le nombre de clients

	cout = matrix_cout(data) # le cout de transport du client i (0:n) au j (0:n) (l'indice 1 est le centre de depot) 	

	# Création d'un modèle. Ce modèle fera l'interface avec le solveur GLPK
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
			@constraint(m, (demande[i,t]-(Q + demande[i,t])*(1 - x[i,j]))<= (w[i]-w[j])) # Contrainte (10)
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

		#Résolution du problème d'optimisation linéaire m par le solveur GLPK
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

# Tests
# pathFileData = "PRP_instances/B_200_instance30.prp"

# Récupérer le q dans LSP
# p, y, I, q = PL_LSP(pathFileData, 0, false)

# x, w = PL_VRP(pathFileData, q, 2, false) # Le q vient de LSP.jl
# println("x=", x)
# println("w=", w)





