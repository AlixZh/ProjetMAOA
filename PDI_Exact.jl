using JuMP
using CPLEX

include("PDI2.jl")
include("PDI1.jl")
include("toolsPDI.jl")
include("tools.jl")
include("BC.jl")




function PDI_exact(data,numAlgo)
	"""
	renvoie la solution exacte sans Branch and cut
	trouve la solution que pour des instances data de taille petite
	renvoie true si il y a une solution optimale
	"""
	if(numAlgo==1)
		m=PLNE_PDI1(data)
	elseif(numAlgo==2)
		m=PLNE_PDI2(data)
	end
	#Résolution du problème d'optimisation linéaire m par le solveur GLPK
	println("Résolution par le solveur linéaire choisi")
	optimize!(m)
	println()

	# Mais on peut vouloir récupérer juste une information précise
	println("Récupération et affichage \"à la main\" d'informations précises")
	status = termination_status(m)
	if status == JuMP.MathOptInterface.INFEASIBLE
		println("Le problème n'est pas réalisable")
	elseif status == JuMP.MathOptInterface.DUAL_INFEASIBLE
		println("Le problème est non borné")
	elseif status == JuMP.MathOptInterface.OPTIMAL
		# Affiche tous les détails d'une solution à l'écran
		#println("Affichage de tous les détails de la solution avec la commande solution_summary")
		#println(solution_summary(m, verbose=true))
		println()
		println("Solution du PLNE PDI ",numAlgo)
		println("Valeur optimale = ", objective_value(m))
		#println("Solution primale optimale :")
		#println("\t x = ", value.(x))
		println("Temps de résolution :", solve_time(m))
		if(numAlgo==1)
			xopt,yopt,qopt=get_param_opt1(data,m)
		elseif(numAlgo==2)
			xopt,yopt,qopt=get_param_opt2(data,m)
		end
		
		return true,objective_value(m), solve_time(m) ,xopt,yopt,qopt
	else
		println("Problème lors de la résolution")
	end
	return false,0,0,0,0,0
	
end


