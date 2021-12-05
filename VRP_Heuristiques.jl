# using JuMP
# using CPLEX
using Base.Iterators
using LinearAlgebra

include("tools.jl") # Pour extraire les données des fichiers
include("LSP.jl") # Pour pouvoir utiliser le q

# const OPTIMAL = JuMP.MathOptInterface.OPTIMAL
# const INFEASIBLE = JuMP.MathOptInterface.INFEASIBLE
# const UNBOUNDED = JuMP.MathOptInterface.DUAL_INFEASIBLE;

function Bin_Packing(data, demande, t)
	"""
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	demande un tableau [1:n, 1:l] de taille n*l, demande[i,t] est la demande du client i au temps t, soit q[i,t] dans la notation LSP
	t l'instant de temps à considérer
	"""

	# Récupération des données
    Q = data["Q"] # Q la capacité maximale de chaque véhicule
	n = data["n"] # n le nombre de clients donc le nombre de noeuds du graphe (sans compter le dépot) ici

	somme = 0 # Poids cumulé dans la boite courante
	boites = Vector{Int64}[] # Ensemble des boites

	#boiteCourante = Vector{Int64}() #PB
	boiteCourante = [0] # Toutes les boites doivent contenir le centre de dépot 0

	for i in 1:n 
		somme += demande[i,t]

		if somme > Q
			#push!(boiteCourante, 0) PB pq ?
			# On enregistre la boite avant d'ajouter l'objet qui fait déborder la capacité
			push!(boites, boiteCourante)

			# On met l'objet qui a fait déborder la boite précédente dans la nouvelle boite 
			# On met donc la somme courante (du poids dans la boite courante) à jour
			somme = demande[i,t]
			#boiteCourante = Vector{Int64}() PB
			boiteCourante = [0]
		end
		push!(boiteCourante, i)
	end

	# pas besoin ???? PB
	#push!(boiteCourante, 0)
	if length(boiteCourante) > 1
		push!(boites, boiteCourante)
	end

	return boites
end

function Clark_Wright(data, demande, t)
	"""
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	demande un tableau [1:n, 1:l] de taille n*l, demande[i,t] est la demande du client i au temps t, soit q[i,t] dans la notation LSP
	t l'instant de temps à considérer
	"""
	# PB  Il est aussi possible d’utiliser sij = c0i + c0j − αcij + β|c0i + c0j | + γ di+dj / d
	# avec α, β, γ pris entre 0 et 2 afin de prendre un peu en compte le fait
	# que la fusion d´epend plus ou moins des distances et plus ou moins des valeurs des demandes.


	# Récupération des données
    Q = data["Q"] # Q la capacité maximale de chaque véhicule
	n = data["n"] # n le nombre de clients donc le nombre de noeuds du graphe (sans compter le dépot) ici

	cout = matrix_cout(data) # le cout de transport du client i (1:n+1) au j (1:n+1) (l'indice 1 est le centre de depot) 	
	# PB cout indice +1

	# Partons d’une solution à n circuits définis comme les n allers-retours de 0 à chaque client i
	solution = Array[]
	for i in 1:n
		push!(solution, [0,i]) # Un circuit aller-retrour de 0 à i est simplement caractérisé par une boite contenant les deux sommets 0 et i
	end

	# Calcul des s[i,j] 
	s = []
	for i in 1:n
		for j in 1:n
			if i!=j
				# On enregistre dans s à la fois le s[i,j] mais aussi les indices (i,j) pour ne pas les perdre lors du tri
				push!(s, [(i,j), cout[1,i+1] + cout[1,j+1] + cout[i+1,j+1]])
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
		boite_j = [] # PB NaN
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
			demande_fusion=0 
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

#PB sectorielle a revoir
function Sectorielle(data, demande, t, angle) #ON SUPPOSE QUE 360 EST DIVISIBLE PAR ANGLE (sinon trop compliqué, flemme)
	"""
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	demande un tableau [1:n, 1:l] de taille n*l, demande[i,t] est la demande du client i au temps t, soit q[i,t] dans la notation LSP
	t l'instant de temps à considérer
	angle PB
	"""

	function vecteurEntreDeuxPoints(pt1,pt2)
		return [pt2[1]-pt1[1],pt2[2]-pt1[2]]
	end

	# Récupération des données
    Q = data["Q"] # Q la capacité maximale de chaque véhicule
	coord = data["coord"] # coord les coordonnees des positions des clients et du centre de depot

	cout = matrix_cout(data) # le cout de transport du client i (1:n+1) au j (1:n+1) (l'indice 1 est le centre de depot) 	
	# PB cout indice +1 

	distance_du_point_le_plus_eloigne=0
	origin = coord[0]
	#Calcul de la distance du point le plus eloigné
	for node in coord
		distance=sqrt((node["x"]-origin["x"])^2+(node["y"]-origin["y"])^2)
		if(distance_du_point_le_plus_eloigne<distance)
			distance_du_point_le_plus_eloigne=distance
		end
	end

	
	#Calcul distance point fictif a l'origine
		#La distance de ces points fictif à l'origine ne peut pas être égale 
		#a la distance du point le plus eloigné car celui ci peut ne pas être aucun des triangles
		#C'est pourquoi que la distance des points fictifs à l'origine sera égal à
		#la somme de la distance du point le plus eloigné + la distance entre le (point D qui est le point au milieu
		# du segment entre deux points fictifs quelconquz) et le
		#(point E de fin du segment qui commence à l'origine, passant par D ) 

	distance_point_fictif_a_lorigine=distance_du_point_le_plus_eloigne*(2-cos(angle/2))

	#Création des points fictifs pour créer les secteurs en triangles
	angleNormalise=angle/360
	triangles=[]
	for i in 1:360/angle
		push!(triangles,[(float(origin["x"]),float(origin["y"]))])
	end
	k=1
	for i in 1:2*360/angle
		point=(origin["x"]+(1-(i-1)*angleNormalise)*distance_point_fictif_a_lorigine,origin["y"]+(i-1)*angleNormalise*distance_point_fictif_a_lorigine)
		push!(triangles[k],point)
		k+=Integer((i+1)%2)
	end

	#Création secteurs (il y en a 360/angle)
	secteurs=[]
	for i in 1:length(triangles)
		push!(secteurs,[origin])
	end
	initial=true
	for node in coord
		if initial
			initial=false
			continue
		end
		k=1
		for (Ori,A,B) in triangles #detection de quel secteur est situé le point
			# AB=vecteurEntreDeuxPoints(A,B)
			# BA=vecteurEntreDeuxPoints(B,A)

			# OriA=vecteurEntreDeuxPoints(Ori,A)
			# AOri=vecteurEntreDeuxPoints(A,Ori)
			
			# OriB=vecteurEntreDeuxPoints(Ori,B)
			# BOri=vecteurEntreDeuxPoints(B,Ori)
			
			# ANode=vecteurEntreDeuxPoints(A,(node["x"],node["y"]))
			# BNode=vecteurEntreDeuxPoints(B,(node["x"],node["y"]))
			# OriNode=vecteurEntreDeuxPoints(C,(node["x"],node["y"]))
			
			# #print(AB,BA,OriA,AOri,OriB,BOri,ANode,BNode,OriNode)
			# if dot(cross(AB,ANode),cross(ANode,AOri))>=0 && dot(cross(BA,BNode),cross(BNode,BOri))>=0 && dot(cross(OriA,OriNode),cross(OriNode,OriB))>=0
			# 	push!(secteurs[k],node)
			# 	break
			# end

			cond1=(A[1]-node["x"])*(B[2]-node["y"])-(A[2]-node["y"])*(B[1]-node["x"])
			cond2=(B[1]-node["x"])*(Ori[2]-node["y"])-(B[2]-node["y"])*(Ori[1]-node["x"])
			cond3=(Ori[1]-node["x"])*(A[2]-node["y"])-(Ori[2]-node["y"])*(A[1]-node["x"])
			if((cond1>=0 && cond2>=0 && cond3>=0) || (cond1<0 && cond2<0 && cond3<0))
				push!(secteurs[k],node)
				break
			end
			k+=1
		end
	end

	ens_de_ens_de_circuits=[]
	for secteur in secteurs
		push!(ens_de_ens_de_circuits,Clark_Wright(data,secteur,demande,t)) #PB
	end

	ens_final=[]
	for ens_de_circuits in ens_de_ens_de_circuits
		for circuit in ens_de_circuits
			push!(ens_final,circuit)
		end
	end
	return ens_final
end

function Boites_heuristique(data, demande, t, heuristique)
	"""
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	demande un tableau [1:n, 1:l] de taille n*l, demande[i,t] est la demande du client i au temps t, soit q[i,t] dans la notation LSP
	t l'instant de temps à considérer
	heuristique une chaine de caractère déterminant quelle heuristique on va utiliser
	"""
	if heuristique == "BP"
		return Bin_Packing(data, demande, t)
	end

	if heuristique == "CW"
		return Clark_Wright(data, demande, t)
	end

	# if heuristique == "Sec" #PB
	# 	return Sectorielle(data, demande, t)
	# end
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
	# PB cout indice +1

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

function Cout_Circuit(data, circuit)
	"""
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	circuit un circuit représenté sous forme de liste 
	"""

	# Récupération des données
    cout = matrix_cout(data) # le cout de transport du client i (1:n+1) au j (1:n+1) (l'indice 1 est le centre de depot) 	
	# PB cout indice +1

	cout_total = 0

	for (from, to) in zip(circuit[1:end-1], circuit[2:end])
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
	# PB cout indice +1

	cout_total = 0

	for circuit in circuits
		# On ajoute le coût du circuit que l'on traite
		cout_total += Cout_Circuit(data, circuit)	
	end

	return cout_total
end

# PB méthode itérative (métaheuristique)
function VRP_iteratif(data, demande, t, heuristique)
	"""
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	demande un tableau [1:n, 1:l] de taille n*l, demande[i,t] est la demande du client i au temps t, soit q[i,t] dans la notation LSP
	t l'instant de temps à considérer
	heuristique une chaine de caractère déterminant quelle heuristique on va utiliser
	"""
	# Récupération des données
    k = data["k"] # k le nombre de véhicules
	n = data["n"] # n le nombre de clients

	cout = matrix_cout(data) # le cout de transport du client i (1:n+1) au j (1:n+1) (l'indice 1 est le centre de depot) 	
	# PB cout indice +1

	# A l’issue d’une étape gloutonne, on obtient une solution réalisable ou alors une solution utilisant plus de k tournées (chacune réalisable)
	boites = Boites_heuristique(data, demande, t, heuristique)

	if length(boites) < k
		# On ajoute des tournées vides pour permettre des transition vers ces nouvelles tournées
		len = length(boites)
		for i in 1:k-len
			push!(boites, [0])
		end
	end

	if length(boites) > k
		# PB
		
		# Si le nombre de tourn´ees de la solution initiale est sup´erieur `a m, il est possible d’utiliser
		# dans un premier temps une fonction objective artificielle pour forcer `a r´eduire le nombre de tourn´ees.
		
	end

	# A partir d'ici, le nombre de tournées est exactement égal à k
	nbIteMax = 10

	nbIte = 0
	changement = true # Vaut true si on a changé quelque chose lors de l'itération courrante
    while nbIte <= nbIteMax && changement
		changement = false

        # Possibilité pour un client de changer de tournées (notamment si il y a des tournées vides)
		for client in 1:n
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

			cout_meilleure_boite = Cout_heur_Boite(data, boite_client)
			indice_meilleure_boite = indice_boite_client

			for indice_boite in 1:length(boites) 
				if indice_boite != indice_boite_client
					nouvelle_boite = copy(boites[indice_boite])
					push!(nouvelle_boite, client)
					cout_nouvelle_boite = Cout_heur_Boite(data, nouvelle_boite)
					if cout_meilleure_boite > cout_nouvelle_boite
						cout_meilleure_boite = cout_nouvelle_boite
						indice_meilleure_boite = indice_boite
						changement = true 
					end
				end
			end
			
			# On le change de boite si c'est moins couteux qu'il soit dans celle-ci plutôt que dans sa boîte actuelle
			if indice_meilleure_boite != indice_boite_client

				# Création des deux boites qui ont changé lorsque le client change de boite (son ancienne et sa nouvelle boite)
				nouv_meilleure_boite = copy(boites[indice_meilleure_boite])
				push!(nouv_meilleure_boite, client)
				nouv_boite_sans_client = []
				for elem in boite_client
					if elem != client
						push!(nouv_boite_sans_client, elem)
					end
				end

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
        nbIte += 1
    end

	#PB Une m´etaheuristique it´erative peut alors ˆetre utilis´ee sur la base de plusieurs voisinages:
	# - supprimer une tourn´ee vide

	# PB il faut ensuite faire le LSP pour résoudre VRP heuristiquement 
	# (les heuristiques ci-dessus partitionnent juste les clients)
	# Permet, grâce à la résolution du problème de voyageur de commerce (TSP), de passer de la partition des clients aux tournées
	# On veut optimiser chacune des tournées (càd remettre dans un ordre 'optimal' chacune des boites)

	for boite in boites # Pour chaque boite on regarde les voisinages classiques du TSP (2-opt) pour améliorer chaque tournée indépendemment
		if length(boite)>3 # Sinon, on ne peut pas faire de voisinage 2-opt
			for ind_elem1 in 1:length(boite)
				for ind_elem2 in elem1+1:length(boite)
					if ind_elem1 != ind_elem2 && ind_elem1+1 != ind_elem2 && ind_elem1-1 != ind_elem2
						if ind_elem1 == 1 && ind_elem2 == length(boite)
							continue
						end

					end
				end
			end
		end
	end

	return boites
end

# Tests
pathFileData = "PRP_instances/B_200_instance30.prp"
data = Read_file(pathFileData)
t = 2

# Récupérer le q dans LSP
p, y, I, q = PL_LSP(data, 0, false)

# Heuristique Bin Packing
boites_BP = Bin_Packing(data, q, t) # Le q vient de LSP.jl
println("Résultat avec Bin Packing :")
println(boites_BP)
println("Nombre de tournées = ", length(boites_BP)) 

# Heuristique Clark-Wright
boites_CW = Clark_Wright(data, q, t) # Le q vient de LSP.jl
println("Résultat avec Clarck Wright :")
println(boites_CW)
println("Nombre de tournées = ", length(boites_CW)) 

# Heuristique sectorielle #PB
# boites_Sec = Sectorielle(data, q, t) # Le q vient de LSP.jl
# println("Résultat avec l'heuristique sectorielle :")
# println(boites_Sec)
# println("Nombre de tournées = ", length(boites_Sec)) 

# VRP heuristique
boites = Boites_heuristique(data, q, t,"BP")
println(boites==boites_BP)

boites = Boites_heuristique(data, q, t,"CW")
println(boites==boites_CW)

# boites = Boites_heuristique(data, q, t,"Sec") #PB
# println(boites==boites_Sec)

# Le coût des circuits obtenus
println("Coût des circuits : ", Cout_Circuits(data, boites))

# Passer des boîtes à des circuits grâce au TSP-> plutot dans itératif ?
#Boites_to_VRP(data, q, t, "CW")

#PB il faut ensuite faire méthode itérative
boites = VRP_iteratif(data, q, t, "CW")
println("Coût des circuits : ", Cout_Circuits(data, boites))

# PB changer pathfile pour ne pas avoir a relire plusieurs fois le fichier