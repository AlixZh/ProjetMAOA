# using JuMP
# using CPLEX
using Base.Iterators
using LinearAlgebra

#include("tools.jl") # Pour extraire les données des fichiers
#include("LSP.jl") # Pour pouvoir utiliser le q
include("TSP.jl")

# const OPTIMAL = JuMP.MathOptInterface.OPTIMAL
# const INFEASIBLE = JuMP.MathOptInterface.INFEASIBLE
# const UNBOUNDED = JuMP.MathOptInterface.DUAL_INFEASIBLE;

function Bin_Packing(data, demande, t)
	"""
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	demande un tableau [1:n, 1:l] de taille n*l, demande[i,t] est la quantité produite pour le revendeur i à la période t, soit q[i,t] dans la notation LSP
	t l'instant de temps à considérer
	"""

	# Récupération des données
    Q = data["Q"] # Q la capacité maximale de chaque véhicule
	n = data["n"] # n le nombre de clients donc le nombre de noeuds du graphe (sans compter le dépot) ici

	somme = 0 # Poids cumulé dans la boite courante
	boites = Vector{Int64}[] # Ensemble des boites

	boiteCourante = [0] # Toutes les boites doivent contenir le centre de dépot 0

	for i in 1:n 
		if demande[i,t] != 0 # On ne traite que les clients qui ont une demande pour ce pas de temps
			somme += demande[i,t]

			if somme > Q
				# On enregistre la boite avant d'ajouter l'objet qui fait déborder la capacité
				push!(boites, boiteCourante)

				# On met l'objet qui a fait déborder la boite précédente dans la nouvelle boite 
				# On met donc la somme courante (du poids dans la boite courante) à jour
				somme = demande[i,t]
				boiteCourante = [0]
			end
			push!(boiteCourante, i)
		end
	end

	if length(boiteCourante) > 1
		push!(boites, boiteCourante)
	end

	return boites
end

function Clark_Wright(data, demande, t)
	"""
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	demande un tableau [1:n, 1:l] de taille n*l, demande[i,t] est la quantité produite pour le revendeur i à la période t, soit q[i,t] dans la notation LSP
	t l'instant de temps à considérer
	"""
	# PB  Il est aussi possible d’utiliser sij = c0i + c0j − αcij + β|c0i + c0j | + γ di+dj / d
	# avec α, β, γ pris entre 0 et 2 afin de prendre un peu en compte le fait
	# que la fusion d´epend plus ou moins des distances et plus ou moins des valeurs des demandes. A tester
	alpha = 1
	beta = 1
	gamma = 1

	# Récupération des données
    Q = data["Q"] # Q la capacité maximale de chaque véhicule
	n = data["n"] # n le nombre de clients donc le nombre de noeuds du graphe (sans compter le dépot) ici
	cout = matrix_cout(data) # le cout de transport du client i (1:n+1) au j (1:n+1) (l'indice 1 est le centre de depot) 	
	# Pour rappel, tous les indices de cout doivent être augmentés de 1 par rapport aux indices de l'énoncé (en julia, les indices commencent à 1)

	# On créer l'ensemble de client à traiter
	clients_t = []
	d = 0
	for i in 1:n
		if demande[i,t] != 0 # On ne traite que les clients qui ont une demande pour ce pas de temps
			d += demande[i,t]
			push!(clients_t, i)
		end
	end

	# Partons d’une solution à n circuits définis comme les n allers-retours de 0 à chaque client i
	solution = Array[]
	for i in clients_t 
		push!(solution, [0,i]) # Un circuit aller-retrour de 0 à i est simplement caractérisé par une boite contenant les deux sommets 0 et i
	end

	# Calcul des s[i,j] 
	s = []
	for i in clients_t 
		for j in clients_t 
			if i!=j
				# On enregistre dans s à la fois le s[i,j] mais aussi les indices (i,j) pour ne pas les perdre lors du tri
				push!(s, [(i,j), cout[1,i+1] + cout[1,j+1] - alpha*cout[i+1,j+1] + beta*abs(cout[1,i+1] + cout[1,j+1]) + gamma*(demande[i,t]+demande[j,t])/d]) #PB divisé par d ?
			end
		end
	end

	# Trier s par s[i,j] décroissant
	sort!(s, by = x->x[2], rev=true) 

	# Parcours dans l'ordre des s[i,j] décroissant pour la fusion des circuits (boites)
	for k in 1:length(s)
		# On récupère les sommets i et j associés au s[i,j] que l'on traite
		i = s[k][1][1]
		j = s[k][1][2]
		boite_i = []
		boite_j = [] 
		indice_boite_i = -1
		indice_boite_j = -1

		# On cherche les circuits (boites) dans lesquels sont i et j
		for indice_boite in 1:length(solution) 
			if i in solution[indice_boite]
				boite_i = solution[indice_boite]
				indice_boite_i = indice_boite
			end

			if j in solution[indice_boite]
				boite_j = solution[indice_boite]
				indice_boite_j = indice_boite
			end

			if indice_boite_i != -1 && indice_boite_j != -1 # On a donc trouvé les deux boites que l'on cherchait
				break
			end
		end

		# Si i et j ne sont pas dans la même boite
		if indice_boite_j != indice_boite_i 
			# On crée la boite constituée des éléments de la boite de i auquel on ajoute tous les éléments de la boite de j (fusion des boites)
			boite_fusionnee = copy(boite_i) 
			for elem in boite_j
				if !(elem in boite_fusionnee)
					push!(boite_fusionnee, elem)
				end
			end

			# On calcule la demande du nouveau circuit caractérisé par la boite boite_fusionnee
			demande_fusion = 0 
			for l in 2:length(boite_fusionnee) # Commence à 2 car les indices commencent à 1 et que le dépot est le premier élément et n'a pas de demande
				demande_fusion += demande[boite_fusionnee[l], t]
			end

			# Si la réunion des deux circuits ne dépasse pas la demande transportable par un seul véhicule, réunir en un seul circuit les deux circuits
			if demande_fusion<=Q
				# On supprime les circuits associés aux boites boite_i et boite_j
				nouv_solution = []
				for indice_boite in 1:length(solution)
					if indice_boite == indice_boite_i || indice_boite == indice_boite_j 
						continue # On ne garde pas les anciens circuits (boites) dans lesquels étaient i ou j
					end
					push!(nouv_solution, solution[indice_boite])
				end
				solution = nouv_solution
				push!(solution, boite_fusionnee) # On ajoute également le circuit constitué de la fusion des deux boites
			end
		end
	end

	return solution 
end

function Sectorielle(data, demande, t) 
	"""
	Cette heuristique n’a de sens que lorsque le dépôt est plus ou moins au centre du nuage des clients (ce qui est souvent le cas)

    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	demande un tableau [1:n, 1:l] de taille n*l, demande[i,t] est la quantité produite pour le revendeur i à la période t, soit q[i,t] dans la notation LSP
	t l'instant de temps à considérer
	"""

	# Récupération des données
    Q = data["Q"] # Q la capacité maximale de chaque véhicule
	n = data["n"] # n le nombre de clients
	coord = data["coord"] # coord les coordonnees des positions des clients et du centre de depot (commence à 1) 
	# Pour rappel, tous les indices de coord doivent être augmentés de 1 par rapport aux indices de l'énoncé (en julia, les indices commencent à 1)

	centre = [coord[1][1], coord[1][2]]
	repere = [1, 0]

	# On créer l'ensemble de client à traiter
	clients_t = []
	for i in 1:n
		if demande[i,t] != 0 # On ne traite que les clients qui ont une demande pour ce pas de temps
			push!(clients_t, i)
		end
	end

	clients_angle = []
	for i in clients_t # On parcourt tous les noeuds à traiter à ce pas de temps
		# On centre les coordonnées, càd on màj les coordonnées en prennant en compte que le point (0,0) est le centre de dépôt
		point_centre = [coord[i+1][1], coord[i+1][2]] - centre 
		
		# Calcul du degré entre le vecteur directionnel [1, 0] et le vecteur centre vers le point (modulo 360 degrés)
		degre = mod(rad2deg(atan(repere...) - atan(point_centre...)), 360)
		push!(clients_angle, (i, degre))
	end

	# Trier clients_angle par angle avec le vecteur directionnel [1, 0] croissant
	sort!(clients_angle, by = x->x[2]) 

	circuits = Vector{Int64}[]
	circuit = [0]
	somme = 0

	for cli_ang in clients_angle
		cli = cli_ang[1]
		demand = demande[cli, t]
		somme += demand

		if somme > Q
			push!(circuits, circuit)

			# Réinitialisation
			somme = demand
			circuit = [0]
		end

		push!(circuit, cli)
	end

	# On ajoute le dernier circuit à la liste de circuits
	if length(circuit) > 1
		push!(circuits, circuit)
	end

	return circuits
end

function Boites_heuristique(data, demande, t, heuristique)
	"""
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	demande un tableau [1:n, 1:l] de taille n*l, demande[i,t] est la quantité produite pour le revendeur i à la période t, soit q[i,t] dans la notation LSP
	t l'instant de temps à considérer
	heuristique une chaine de caractère déterminant quelle heuristique on va utiliser (ou "PL" si on veut résoudre exactement avec le PL)
	"""

	if heuristique == "PL" 
		x, w = PL_VRP(data, demande, t, false)
		return VRP_to_Circuit(data, x) 
	end
	
	if heuristique == "BP"
		return Bin_Packing(data, demande, t)
	end

	if heuristique == "CW"
		return Clark_Wright(data, demande, t)
	end

	if heuristique == "Sec"
		return Sectorielle(data, demande, t)
	end
end

function Cout_heur_Boite(data, boite)
	"""
	Retourne une borne inférieure sur le coût de la tournée en entrée
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	boite une boite, càd la liste des sommets d'une tournée 
	"""
				
	# Récupération des données
    cout = matrix_cout(data) # le cout de transport du client i (1:n+1) au j (1:n+1) (l'indice 1 est le centre de depot) 	
	# Pour rappel, tous les indices de cout doivent être augmentés de 1 par rapport aux indices de l'énoncé (en julia, les indices commencent à 1)

	cout_total = 0

	for i in boite
		coutMin1 = maximum(cout[:])
		coutMin2 = maximum(cout[:])
		for j in boite
			if i != j
				if cout[i+1, j+1]<coutMin1
					if coutMin2 > coutMin1
						coutMin2 = coutMin1
					end
					coutMin1 = cout[i+1, j+1]
				end
				if cout[i+1, j+1]>coutMin1 && cout[i+1, j+1]<coutMin2
					coutMin2 = cout[i+1, j+1]
				end
			end
		end
		cout_total += coutMin1 + coutMin2
	end

	return cout_total/2
end

function MinCout_Ajout(cout, boite, i) 
	"""
	Retourne la valeur du coût minimum de la boîte dans laquelle on a insérer i (et l'indice auquel il faut l'insérer pour atteindre ce coût)
    Paramètres : 
    cout la matrice de cout de transport, cout[i,j] du client i (1:n+1) au j (1:n+1) (l'indice 1 est le centre de depot)
	boite une boite, càd la liste des sommets d'une tournée 
	i un élément (numéro de noeud) à insérer
	"""

	# Initialisation du coût
	circuit = copy(boite)
	push!(circuit, i)
	cout_min = Cout_Circuit(cout, circuit)
	indice_min = length(boite)+1

	for indice in 2:length(boite) # Commence à 2 car l'indice 1 est réservé au centre de dépôt
		circuit = copy(boite)
		insert!(circuit, indice, i)  
		cout_courrant = Cout_Circuit(cout, circuit)

		if cout_min > cout_courrant
			cout_min = cout_courrant
			indice_min = indice
		end
	end

	return cout_min, indice_min 
end

function Cout_Circuit(cout, circuit)
	"""
    Paramètres : 
    cout la matrice de cout de transport, cout[i,j] du client i (1:n+1) au j (1:n+1) (l'indice 1 est le centre de depot)
	circuit un circuit représenté sous forme de liste 
	"""
	cout_total = 0

	# println("PB n=",length(cout))
	# println("PB circuit :", circuit)
	for (from, to) in zip(circuit[1:end-1], circuit[2:end])
		# println("PB from : ", from)
		# println("PB to : ", to)
		# println("PB vecteur cout", cout)
		# println("PB cout ", cout[from+1, to+1])
		cout_total += cout[from+1, to+1]
	end

	cout_total += cout[circuit[end]+1, circuit[1]+1] # Sans oublier le coût de l'arc pour revenir au dépôt
	
	return cout_total
end

function Cout_Circuits(data, circuits)
	"""
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	circuits un ensemble de circuits représenté sous forme de liste de liste
	"""

	# Récupération des données
    cout = matrix_cout(data) # le cout de transport du client i (1:n+1) au j (1:n+1) (l'indice 1 est le centre de depot) 	
	# Pour rappel, tous les indices de cout doivent être augmentés de 1 par rapport aux indices de l'énoncé (en julia, les indices commencent à 1)

	cout_total = 0

	for circuit in circuits
		#println("PB circuit ",circuit," de cout ", Cout_Circuit(cout, circuit))
		# On ajoute le coût du circuit que l'on traite
		cout_total += Cout_Circuit(cout, circuit)	
	end

	return cout_total
end
 
function VRP_iteratif(data, demande, t, heuristique)
	"""
	Méthode itérative (métaheuristique)
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	demande un tableau [1:n, 1:l] de taille n*l, demande[i,t] est la quantité produite pour le revendeur i à la période t, soit q[i,t] dans la notation LSP
	t l'instant de temps à considérer
	heuristique une chaine de caractère déterminant quelle heuristique on va utiliser (ou "PL" si on veut résoudre exactement avec le PL)
	"""
	# Récupération des données
    k = data["k"] # k le nombre de véhicules
	n = data["n"] # n le nombre de clients
    Q = data["Q"] # Q la capacité maximale de chaque véhicule

	cout = matrix_cout(data) # le cout de transport du client i (1:n+1) au j (1:n+1) (l'indice 1 est le centre de depot) 	
	# Pour rappel, tous les indices de cout doivent être augmentés de 1 par rapport aux indices de l'énoncé (en julia, les indices commencent à 1)

	# On créer l'ensemble de client à traiter
	clients_t = []
	for i in 1:n
		if demande[i,t] != 0 # On ne traite que les clients qui ont une demande pour ce pas de temps
			push!(clients_t, i)
		end
	end

	# A l’issue d’une étape gloutonne, on obtient une solution réalisable ou alors une solution utilisant plus de k tournées (chacune réalisable)
	boites = Boites_heuristique(data, demande, t, heuristique)
	# PB tests -------------------------------------
	# function histogram(s)
	# 	d = Dict()
	# 	for c in s
	# 		if !(c in keys(d))
	# 			d[c] = 1
	# 		else
	# 			d[c] += 1
	# 		end
	# 	end
	# 	return d
	# end

	# pb = false
	# for boite in boites
	# 	if boite[1] != 0
	# 		pb = true
	# 	end
	# 	dict = histogram(boite)
	# 	for k in keys(dict)
	# 		if dict[k]>1
	# 			println("PB boite il y a ", dict[k], " numéros ", k)
	# 			pb = true
	# 		end
	# 	end
	# end 
	#println("PB boites = ", pb)
	# PB tests ------------------------------------- 

	# PB modifier dans overleaf que on fait LSP à ce moment là et pas 2fois ou quoi

	# Il faut ensuite faire le TSP pour résoudre le VRP heuristiquement (les heuristiques ci-dessus partitionnent juste les clients)
	# Cela permet, grâce à la résolution du problème de voyageur de commerce (TSP), de passer de la partition des clients aux tournées
	# On veut optimiser chacune des tournées (càd remettre dans un ordre 'optimal' chacune des boites)

	# PB marche pas :
	# On résout le problème du voyageur de commerce (TSP) limité aux sommets de chaque boîte par PLNE
	# new_boites = []
	# for boite in boites
	# 	println("-----------------------------------------------------------------")
	# 	println("PB boite avant TSP = ", boite)
	# 	x = PL_TSP(data, boite, false)
	# 	new_boite = TSP_to_Circuit(data, boite, x)
	# 	push!(new_boites, new_boite) 
	# 	println("PB boite après TSP = ", new_boite)
	# 	println("-----------------------------------------------------------------")
	# end
	# boites = new_boites 

	# PB(par une heuristique gloutonne au plus proche sommet en démarrant du dépôt)
	function TSP_heur(boites, co)
		new_boites = []
		for boite in boites
			new_boite = [0]
			filter!(e->e!=0, boite)
			while length(boite) > 1
				# On cherche le sommet le plus proche du dernier élément parmi ceux qu'il reste
				CoutMin = co[new_boite[end]+1, boite[1]+1]
				elemMin = boite[1]
				for elem in boite
					if CoutMin > co[new_boite[end]+1, elem+1]
						CoutMin = co[new_boite[end]+1, elem+1]
						elemMin = elem
					end
				end
				push!(new_boite, elemMin)
				filter!(e->e!=elemMin, boite)
			end
			if length(boite) == 1
				push!(new_boite, boite[1]) # On ajoute le dernier élément
			end
			push!(new_boites, new_boite)
		end
		return new_boites
	end

	boites = TSP_heur(boites, cout)

	# PB tests -------------------------------------
	# pb = false
	# for boite in boites
	# 	if boite[1] != 0
	# 		pb = true
	# 	end
	# # 	dict = histogram(boite)
	# # 	for k in keys(dict)
	# # 		if dict[k]>1
	# # 			println("PB boite il y a ", dict[k], " numéros ", k)
	# # 			pb = true
	# # 		end
	# # 	end
	# end 
	#println("PB boites apres TSP glouton = ", pb)
	# PB tests -------------------------------------

	# println("PB Coût après TSP : ", Cout_Circuits(data, boites)) 


	#println("PB AVANT BOITE, cout = ", Cout_Circuits(data, boites))
	# println("PB BOITES = ", boites)
	# if length(boites) < k # PB a decommenter
	# 	# On ajoute des tournées vides pour permettre des transition vers ces nouvelles tournées
	# 	len = length(boites)
	# 	for i in 1:k-len
	# 		push!(boites, [0])
	# 	end
	# end

	if length(boites) < k #PB ajouter dans le overleaf cette initialisation
		# On sépare les tournées jusqu'à utiliser tous les véhicules
		len = length(boites)
		for i in 1:k-len # On ajoute k-len boîtes pour atteindre k boîtes
			# On choisit la boîte comportant le plus d'éléments
			ind_grande_boite = 1
			for ind_boite in 2:length(boites)
				if length(boites[ind_boite]) > length(boites[ind_grande_boite])
					ind_grande_boite = ind_boite
				end
			end

			# On sépare la boîte en deux : notons que la capacité ne dépasse forcément pas celle du camion puisque ça ne dépassait pas 
			# pour la boîte entière et que l'heuristique TSP dans laquelle on prend le sommet le moinsloin à chaque fois est vérifiée 
			# pour ces sous ensembles, nous n'avons donc pas besoin de le refaire
			boite_actu = boites[ind_grande_boite]
			tab1 = boite_actu[1:floor(Int,length(boite_actu)/2)] # Commence bien par le dépôt 0
			tab2 = boite_actu[floor(Int,length(boite_actu)/2)+1:length(boite_actu)]
			boites[ind_grande_boite] = tab1
			# Il faut que tab2 commence par le centre de dépôt
			insert!(tab2, 1, 0)
			push!(boites, tab2)
		end
	end

	#println("PB BOITES APRES Ajout boites vides = ", boites)
	#println("PB Coût après ajout de boites : ", Cout_Circuits(data, boites))

	if length(boites) > k # PB ajouter dans le overleaf cette initialisation
		# On réunit les tournées jusqu'à n'utiliser pas plus de k véhicules
		len = length(boites)
		for i in 1:len-k # On enlève len-k boîtes pour atteindre k boîtes
			changementNonEffectue = true
			nbIteEchec = 0
			while changementNonEffectue
				if nbIteEchec == 0
					# On choisit les deux boîtes comportant le moins d'éléments
					ind_petite_boite = 1
					ind_petite_boite2 = 1
					for ind_boite in 2:length(boites)
						if length(boites[ind_boite]) < length(boites[ind_petite_boite])
							if length(boites[ind_petite_boite2]) > length(boites[ind_petite_boite])
								ind_petite_boite2 = ind_petite_boite
							end
							ind_petite_boite = ind_boite
						end

						if length(boites[ind_boite]) > length(boites[ind_petite_boite]) && length(boites[ind_boite]) < length(boites[ind_petite_boite2])
							ind_petite_boite2 = ind_boite
						end
					end

					# Si la demande des deux boîtes n'excède pas la capacité du camion, c'est bon
					# Il faut également que la capacité du camion (Q) ne soit pas dépassée dans la nouvelle boîte
					nouvelle_demande = 0
					for elemq in boites[ind_petite_boite]
						if elemq != 0
							nouvelle_demande += demande[elemq, t]
						end
					end

					for elemq in boites[ind_petite_boite2]
						if elemq != 0
							nouvelle_demande += demande[elemq, t]
						end
					end
				
					if nouvelle_demande <= Q && ind_petite_boite != ind_petite_boite2
						# On fait TSP sur la nouvelle boîte réunissant les 2 boîtes précédentes
						boite2 = boites[ind_petite_boite2]
						nouv_boite = [boites[ind_petite_boite] ; boite2[2:length(boite2)]] # On ne prend pas deux fois le dépôt

						nouv_boite_tsp = TSP_heur([nouv_boite], cout)[1]
						boites[ind_petite_boite] =  nouv_boite_tsp

						# On retire la seconde boîte
						deleteat!(boites, ind_petite_boite2)
						changementNonEffectue = false
					end
					nbIteEchec += 1 # Ne restera à 1 que si on reboucle
				else
					for ind_petite_boite in 1:length(boites)-1
						if changementNonEffectue == false
							break
						end
						for ind_petite_boite2 in ind_petite_boite+1:length(boites)
							if changementNonEffectue == false
								break
							end
							# Si la demande des deux boîtes n'excède pas la capacité du camion, c'est bon
							# Il faut également que la capacité du camion (Q) ne soit pas dépassée dans la nouvelle boîte
							nouvelle_demande = 0
							for elemq in boites[ind_petite_boite]
								if elemq != 0
									nouvelle_demande += demande[elemq, t]
								end
							end

							for elemq in boites[ind_petite_boite2]
								if elemq != 0
									nouvelle_demande += demande[elemq, t]
								end
							end
						
							if nouvelle_demande <= Q
								# On fait TSP sur la nouvelle boîte réunissant les 2 boîtes précédentes
								boite2 = boites[ind_petite_boite2]
								nouv_boite = [boites[ind_petite_boite] ; boite2[2:length(boite2)]] # On ne prend pas deux fois le dépôt

								nouv_boite_tsp = TSP_heur([nouv_boite], cout)[1]
								boites[ind_petite_boite] =  nouv_boite_tsp

								# On retire la seconde boîte
								deleteat!(boites, ind_petite_boite2)
								changementNonEffectue = false
							end

							if ind_petite_boite == length(boites)-1 && ind_petite_boite2 == length(boites) && changementNonEffectue
								# Dans ce cas aucun changement ne peut être effectué avec un tel voisinnage sans dépasser la taille du camion
								# Dans ce cas on regarde déjà si le problème est réalisable
								if sum(demande[i,t] for i in 1:n) > Q*k # Normalement impossible grâce au LSP
									println("!! Attention problème NON REALISABLE !!")
									return 0
								end

								# Et si il est bien réalisable, on utilise le PL du VRP
								x, w = PL_VRP(data, demande, t, false)
								return VRP_to_Circuit(data, demande, x) # On retourne directement puisque le PL résout exactement le problème

								# Autres possibilités :
								# Si le nombre de tournées de la solution initiale est supérieur à k, il est possible d’utiliser
								# dans un premier temps une fonction objective artificielle pour forcer à réduire le nombre de tournées.

								# Ou on prend les boîtes (en nombre nécessaire pour atteindre k boîtes) de demande la plus faible
								# Puis on répartit dans les autres boîtes jusqu'à dépasser le poids (problème du sac à dos) 
								# mais il n'est pas forcé que l'on trouve une solution réalisable non plus
							end
						end
					end
				end
			end
		end
	end

	# PB tests -------------------------------------
	# pb = false
	# for boite in boites
	# 	if boite[1] != 0
	# 		pb = true
	# 	end
	# # 	dict = histogram(boite)
	# # 	for k in keys(dict)
	# # 		if dict[k]>1
	# # 			println("PB boite il y a ", dict[k], " numéros ", k)
	# # 			pb = true
	# # 		end
	# # 	end
	# end 
	#println("PB boites après ajout circuits vides= ", pb)
	# PB tests -------------------------------------

	# A partir d'ici, le nombre de tournées est exactement égal à k 
	nbIteMax = 5

	nbIte = 0
	changement = true # Vaut true si on a changé quelque chose lors de l'itération courrante
    while nbIte < nbIteMax && changement
		changement = false

        # Possibilité pour un client de changer de tournées (notamment si il y a des tournées vides)
		for client in clients_t
			boite_client = []
			indice_boite_client = -1

			# On cherche le circuit (boite) dans lequel est client
			for indice_boite in 1:length(boites) 
				if client in boites[indice_boite]
					boite_client = boites[indice_boite]
					indice_boite_client = indice_boite
					break
				end
			end

			# PB On prend ici un coût heuristique car on n'a pas encore appliqué un voisinage de TSP dans chacune des boîtes
			# On cherche d'abord à définir les boîtes de coût heuristique les plus faible 
			cout_boite_client = Cout_Circuit(cout, boite_client) # PB Cout_heur_Boite(data, boite_client)

			boite_client_sans = copy(boite_client)
			filter!(e->e!=client, boite_client_sans)
			cout_boite_client_sans = Cout_Circuit(cout, boite_client_sans) # PB Cout_heur_Boite(data, boite_client_sans)

			# println("------PB")
			# println("PB boite_client = ", length(boite_client), " et ", client in boite_client)
			# println("PB boite_client SANS = ", length(boite_client_sans), " et ", client in boite_client_sans)
			cout_meilleur_changement = 0
			indice_meilleure_boite = indice_boite_client
			min_indice_meilleur = 0

			for indice_boite in 1:length(boites) 
				if indice_boite != indice_boite_client

					min_cout, min_indice = MinCout_Ajout(cout, boites[indice_boite], client)
					ancien_cout = cout_boite_client + Cout_Circuit(cout, boites[indice_boite])
					nouveau_cout = cout_boite_client_sans + min_cout

					# Si (le coût de la boîte avec le client + le coût de la boîte d'indice indice_boite 
					# - le coût de l'ancienne boîte du client sans le client - le coût de la boîte d'indice indice_boite avec le client en plus), 
					# càd ce que l'on gagne en changeant le client de boite est plus grand que ce que l'on a trouvé jusqu'alors
					# et que la capacité du camion est toujours repectée, on màj la meilleure boite trouvée
					if ancien_cout - nouveau_cout > cout_meilleur_changement 
						# Il faut également que la capacité du camion (Q) ne soit pas dépassée dans la nouvelle boîte
						nouvelle_demande = 0
						for elemq in boites[indice_boite]
							if elemq != 0
								nouvelle_demande += demande[elemq, t]
							end
						end
						nouvelle_demande += demande[client, t]
					
						if nouvelle_demande <= Q
							cout_meilleur_changement = ancien_cout - nouveau_cout
							indice_meilleure_boite = indice_boite # On garde la boîte qui améliore le plus la solution
							min_indice_meilleur = min_indice # L'indice où il faut insérer client dans boites[indice_meilleure_boite] pour avoir le coût le plus faible
							changement = true 
						end
					end

					# # PB avec heuristique :
					# push!(nouvelle_boite, client)

					# ancien_cout = cout_boite_client + Cout_heur_Boite(data, boites[indice_boite])
					# nouveau_cout = cout_boite_client_sans + Cout_heur_Boite(data, nouvelle_boite)
					
					# # Si (le coût de la boîte avec le client + le coût de la boîte d'indice indice_boite 
					# # - le coût de l'ancienne boîte du client sans le client - le coût de la boîte d'indice indice_boite avec le client en plus), 
					# # càd ce que l'on gagne en changeant le client de boite est plus grand que ce que l'on a trouvé jusqu'alors, on màj la meilleure boite trouvée
					# if ancien_cout - nouveau_cout > cout_meilleur_changement
					# 	cout_meilleur_changement = ancien_cout - nouveau_cout
					# 	indice_meilleure_boite = indice_boite # On garde la boîte qui améliore le plus la solution
					# 	changement = true 
					# end
				end
			end
			
			# On le change de boite si c'est moins couteux qu'il soit dans celle-ci plutôt que dans sa boîte actuelle
			if indice_meilleure_boite != indice_boite_client

				# Création des deux boites qui ont changé lorsque le client change de boite (son ancienne et sa nouvelle boite)
				nouv_meilleure_boite = copy(boites[indice_meilleure_boite])
				insert!(nouv_meilleure_boite, min_indice_meilleur, client)

				nouv_boite_sans_client = copy(boite_client)
				filter!(e->e!=client, nouv_boite_sans_client)

				# On supprime la boite du client et la meilleure boite
				nouv_boites = []
				for indice_boite in 1:length(boites)
					if indice_boite == indice_meilleure_boite  || indice_boite == indice_boite_client
						continue # On ne garde pas les anciennes boites
					end
					push!(nouv_boites, boites[indice_boite])
				end
				boites = nouv_boites
				push!(boites, nouv_meilleure_boite) # On ajoute également les deux boites dans lesquelles client a bougé (était ou sera)
				push!(boites, nouv_boite_sans_client)
			end

		end
		#println("PB Boites =", boites)
        nbIte += 1
		# if changement
		# 	println("PB Coût après mouvement d'un client : ", Cout_Circuits(data, boites)) # PB, "chgmt = ", changement, "heur = ", sum(Cout_heur_Boite(data, boite) for boite in boites))
    	# end
	end

	# PB tests -------------------------------------
	# pb = false
	# for boite in boites
	# 	if boite[1] != 0
	# 		pb = true
	# 	end
	# # 	dict = histogram(boite)
	# # 	for k in keys(dict)
	# # 		if dict[k]>1
	# # 			println("PB boite il y a ", dict[k], " numéros ", k)
	# # 			pb = true
	# # 		end
	# # 	end
	# end 
	#println("PB boites après mouvement client  = ", pb)
	# PB tests -------------------------------------



	
	# PB tests -------------------------------------
	# pb = false
	# for boite in boites
	# 	if boite[1] != 0
	# 		pb = true
	# 	end
	# # 	dict = histogram(boite)
	# # 	for k in keys(dict)
	# # 		if dict[k]>1
	# # 			println("PB boite il y a ", dict[k], " numéros ", k)
	# # 			pb = true
	# # 		end
	# # 	end
	# end 
	#println("PB boites avant mouvement TSP 2-opt = ", pb)
	# PB tests -------------------------------------

	nbIteMaxTSP2opt = 5 # On peut borner le nombre d'itérations à n au carré (ou puissance 4 (n carré arêtes dans un graphe complet) comme introduction de nouvelles arrêtes à chaque itération) où n est le nombre de sommets dans le graphe

	nbIte = 0
	changement = [true for i in 1:length(boites)] # Vaut true si on a changé quelque chose lors de l'itération courrante
    while nbIte < nbIteMaxTSP2opt && maximum(changement) # Si aucun changement sur aucune boîte à l'itération précédente, on sort
		for ind_boite in 1:length(boites) # Pour chaque boite on regarde les voisinages classiques du TSP (2-opt) pour améliorer chaque tournée indépendemment
			if changement[ind_boite] == false # Si aucun changement n'a été effectué sur cette boîte la fois dernière, on passe à la prochaine boîte
				continue
			end

			#println("PB Taille de boite ", ind_boite, " = ", length(boites[ind_boite]))

			changement[ind_boite] = false
			if length(boites[ind_boite])>3 # Sinon, on ne peut pas faire de voisinage 2-opt
				for ind_elem1 in 1:length(boites[ind_boite])-2
					if changement[ind_boite] # On ne fait un changement dans une boîte qu'une fois par itération (sinon problème)
						break
					end
					
					for ind_elem2 in ind_elem1+2:length(boites[ind_boite])
						if changement[ind_boite] # On ne fait un changement dans une boîte qu'une fois par itération (sinon problème)
							break
						end

						if ind_elem1 != ind_elem2 && ind_elem1+1 != ind_elem2 && ind_elem1-1 != ind_elem2
							if ind_elem1 == 1 && ind_elem2 == length(boites[ind_boite]) # Ou l'inverse mais pas possible comme ind_elem2 commence à ind_elem1
								continue
							end

							# Voisin 2-opt de boite en echangeant les arêtes ind_elem1->ind_elem1+1 et ind_elem2->ind_elem2+1 
							# Ces arêtes deviennent ind_elem1->ind_elem2 et ind_elem1+1->ind_elem2+1
							"""ind_elem10 = (ind_elem1+1)%length(boites[ind_boite]) # PB a enlever
							ind_elem12 = (ind_elem1+2)%length(boites[ind_boite])
							ind_elem21 = (ind_elem2+1)%length(boites[ind_boite])
							ind_elem22 = (ind_elem2+2)%length(boites[ind_boite])
							
							# Pour régler le fait que x%x vaut 0 alors que l'on veut que ça vaille x ici
							if ind_elem10 == 0
								ind_elem10 = length(boites[ind_boite])
							end
							if ind_elem12 == 0
								ind_elem12 = length(boites[ind_boite])
							end
							if ind_elem21 == 0
								ind_elem21 = length(boites[ind_boite])
							end
							if ind_elem22 == 0
								ind_elem22 = length(boites[ind_boite])
							end"""
							ind_elem21 = ind_elem2 + 1
							if ind_elem2 == length(boites[ind_boite]) # Pour avoir l'indice suivant du dernier indice
								ind_elem21 = 1
							end

							# println("PB taille boites = ", length(boites))
							# println("PB ind_boite=",ind_boite)
							# println("PB ind_elem1=",ind_elem1)
							# println("PB ind_elem2=",ind_elem2)
							# println("PB boites[ind_boite][ind_elem1] =")
							# println(boites[ind_boite][ind_elem1])
							# println("boites[ind_boite][ind_elem1+1] =")
							# println(boites[ind_boite][ind_elem1+1])
							# println("cout[boites[ind_boite][ind_elem1], boites[ind_boite][ind_elem1+1]] =")
							# println(cout[boites[ind_boite][ind_elem1]+1, boites[ind_boite][ind_elem1+1]+1])
							cout_el10_el12 = cout[boites[ind_boite][ind_elem1]+1, boites[ind_boite][ind_elem1+1]+1]
							cout_el21_el22 = cout[boites[ind_boite][ind_elem2]+1, boites[ind_boite][ind_elem21]+1]
							cout_el10_el21 = cout[boites[ind_boite][ind_elem1]+1, boites[ind_boite][ind_elem2]+1] 
							cout_el12_el22 = cout[boites[ind_boite][ind_elem1+1]+1, boites[ind_boite][ind_elem21]+1] 
								
							if cout_el10_el12 + cout_el21_el22 > cout_el10_el21 + cout_el12_el22  
								#println("PB Taille de boite ", ind_boite, " = ", length(boites[ind_boite]), " avant chgmt")
								# Dans ce cas, c'est une meilleure solution que d'échanger ces deux arcs
								#println("PB CHANGEMENT boite ",ind_boite, " (taille ", length(boites[ind_boite]), ") elem boite[", ind_elem1, "] =", boites[ind_boite][ind_elem1], " et boite[", ind_elem2, "] =", boites[ind_boite][ind_elem2])
													
								#println("PB boite ", ind_boite, " chgmt indices ", ind_elem1, " et ", ind_elem2, " cout: ", cout_el10_el12 + cout_el21_el22, ">", cout_el10_el21 + cout_el12_el22)
								
								
								boite_voisine = copy(boites[ind_boite][1:ind_elem1])
								for indice in ind_elem2:-1:ind_elem1+1
									push!(boite_voisine, boites[ind_boite][indice])
								end

								if ind_elem2 != length(boites[ind_boite]) # Si ind_elem2 est égal à la taille de la liste, on a fini d'écrire boite_voisinne
									for indice in ind_elem21:length(boites[ind_boite])
										push!(boite_voisine, boites[ind_boite][indice])
									end 
								end
								boites[ind_boite] = boite_voisine

								# if boites[ind_boite][1] != 0 # PB s'assurer que 0 au début de la liste
								# 	println("PB Boite ", ind_boite, " ind_elem1 = ", ind_elem1, " ind_elem2 = ", ind_elem2)
								# 	println("PB ", boites[ind_boite])
								# end
								changement[ind_boite] = true # On passe à la prochaine boîte
							end
						end
					end
				end
			end
			# if changement[ind_boite]
			# 	println("PB Coût après TSP 2-opt d'une boite : ", Cout_Circuits(data, boites)) #PB, " modif = ", changement[ind_boite])
			# end
		end
		nbIte += 1
    end

	# PB tests -------------------------------------
	# pb = false
	# for boite in boites
	# 	if boite[1] != 0
	# 		pb = true
	# 	end
	# # 	dict = histogram(boite)
	# # 	for k in keys(dict)
	# # 		if dict[k]>1
	# # 			#println("PB boite il y a ", dict[k], " numéros ", k)
	# # 			pb = true
	# # 		end
	# # 	end
	# end 
	#println("PB boites après mouvement TSP 2-opt = ", pb)
	# PB tests -------------------------------------

	# On retire les circuits vides
	filter!(e->e!=[0], boites)

	return boites
end

# ----- Tests -----

# pathFileData = "PRP_instances/B_200_instance30.prp"
#pathFileData = "PRP_instances/B_100_instance20.prp"
#pathFileData = "PRP_instances/A_100_ABS92_100_5.prp"
#pathFileData = "PRP_instances/A_014_ABS2_15_1.prp"

# data = Read_file(pathFileData)
# t = 2

# Récupérer le q dans LSP
# p, y, I, q = PL_LSP(data, 0, false)

# Heuristique Bin Packing
# boites_BP = Bin_Packing(data, q, t) # Le q vient de LSP.jl
# println("Résultat avec Bin Packing :")
# println(boites_BP)
# println("Nombre de tournées = ", length(boites_BP)) 

# Heuristique Clark-Wright
# boites_CW = Clark_Wright(data, q, t) # Le q vient de LSP.jl
# println("Résultat avec Clark Wright :")
# println(boites_CW)
# println("Nombre de tournées = ", length(boites_CW)) 

# Heuristique sectorielle 
# boites_Sec = Sectorielle(data, q, t) # Le q vient de LSP.jl
# println("Résultat avec l'heuristique sectorielle :")
# println(boites_Sec)
# println("Nombre de tournées = ", length(boites_Sec)) 

# VRP heuristique
# boites = Boites_heuristique(data, q, t,"BP")
# println(boites==boites_BP)

# boites = Boites_heuristique(data, q, t,"CW")
# println(boites==boites_CW)

# boites = Boites_heuristique(data, q, t,"Sec") 
# println(boites==boites_Sec)

# Le coût des circuits obtenus avec l'heuristique CW
# println("Coût des circuits : ", Cout_Circuits(data, boites_CW))

#PB il faut ensuite faire méthode itérative
# boites = VRP_iteratif(data, q, t, "CW") # PB pour nbIteMaxTSP2opt = 10 (9 ok) et pathFileData = "PRP_instances/B_200_instance30.prp"
# println("Coût des circuits : ", Cout_Circuits(data, boites))

# PB traiter le fait qu'on ne prends en compte que les revendeurs à livrer au temps t (dont le q[i,t] != 0), le centraliser
# PB calculer la matrice de cout une fois pour toute (voir Cout_Circuits par exemple, matrix_cout)
# PB utiliser filter!(e->e!="s", liste) plutot que de recopier toute la liste
# PB idem avec insert!(circuit, indice, eleme)  et  insert!(L, indice, elem)
# PB tester et comparer avec les résultats du prof envoyés par mail
# PB bien commenter
# PB dans l'itératif, s'assurer que ça ne dépasse pas la taille du camion Q
# PB permettre d'utiliser tous les camions de façon plus maligne
# PB s'assurer que les circuits commencent à 0 et n'ai aucun doublon

# PB :
# 39218 si on met TSP avant l'ajout des [0], peut-être vaut-il mieux le mettre en haut et considérer les changements avec les vrais couts
# 33091 si on met TSP avant la phase TSP 2-opt
# 26538 pour avant et après
# 21329 ou 24831 avec les capacités respéctées et sans heuristique au début *
# 30157 avec la nouvelle initialisation ou on coupe les boites en 2 plutot que d'ajouter des boites vides 