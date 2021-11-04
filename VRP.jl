using JuMP
using CPLEX


const OPTIMAL = JuMP.MathOptInterface.OPTIMAL
const INFEASIBLE = JuMP.MathOptInterface.INFEASIBLE
const UNBOUNDED = JuMP.MathOptInterface.DUAL_INFEASIBLE;


#le centre/source de dépot a un noeud d'indice 1 
function PL_VRP(lqi,ldi,data)
	"""
	Parametres : 
	lqi : (vecteur), lqi[i] quantité a etre livree au vendeur i, 1<=i<=n
	lqi[i] : le client i doit recevoir lqi[i] 1<i<=n, 
	ldi : les revendeurs qui vont etre livree, ldi[i] : contient le nom du revendeur, 1<=ldi[i]<=n
	1<=i<=length(ldi)
	l'indice 1 est bien le client et non le centre de depot dans ldi
	data : dictionnaire contenant les infos d une instance
	"""

	lci=matrix_cout(data) #l'indice 1 est le centre de depot
						  # le cout du client i et j dans lci : ldi[i+1] et ldi[j+1]
	# Création d'un modèle. Ce modèle fera l'interface avec le solveur GLPK
	m = Model(CPLEX.Optimizer)
	# Création des variables
	@variable(m, 0<=x[1:data["n"]+1 , 1:data["n"]+1]<=1,Int)
	@variable(m,0<=w[1:data["n"]]) #indice 1 correspond au client 1
	@objective(m, Min, sum(x .*lci) )


	@constraint(m,sum(x[1,:])<=data["k"]) #indice 1 correspond au depot 0
	@constraint(m,sum(x[:,1])<=data["k"])
	for i in [1:length(ldi);]
		@constraint(m,sum(x[Int(ldi[i]+1),:])==1)
		@constraint(m,sum(x[:,Int(ldi[i]+1)])==1)
		@constraint(m,0<=w[i]<=data["Q"])
		for j in [1:length(ldi);]
			@constraint(m, (lqi[Int(ldi[i])]-(data["Q"]+lqi[Int(ldi[i])])*(1-x[i+1,j+1])) <= (w[i]-w[j]))
		end
	end

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
	println("Affichage de tous les détails de la solution avec la commande solution_summary")
	println(solution_summary(m, verbose=true))
	println()


	# Mais on peut vouloir récupérer juste une information précise
	println("Récupération et affichage \"à la main\" d'informations précises")
	status = termination_status(m)

	if status == INFEASIBLE
		println("Le problème n'est pas réalisable")
	elseif status == UNBOUNDED
		println("Le problème est non borné")
	elseif status == OPTIMAL
		println("Valeur optimale = ", objective_value(m))
		println("Solution primale optimale :")
		println("\t x = ", value.(x))
		println("Temps de résolution :", solve_time(m))
	else
		println("Problème lors de la résolution")
	end
end




#FONCTIONNE 
dataA_014_ABS1_15_1=Read_file("./PRP_instances/A_014_ABS1_15_1.prp")
qq=dataA_014_ABS1_15_1["d"][:,2]
dd=[1 4 8]
PL_VRP(qq,dd,dataA_014_ABS1_15_1)






