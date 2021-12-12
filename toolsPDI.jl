#---------------------------------------------------------------
#PRP 1
#---------------------------------------------------------------



function PDI1_Exact_to_Circuit(data, x) 
    """
	Permet de représenter la solution renvoyé par le PL sous forme de liste de liste (représentant les circuits, dans l'ordre de passage)
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	x les valeurs des variables binaires xij qui vallent 1 si un véhicule utilise l’arc (i, j) 
	affichage qui vaut true si on veut afficher le graphe
    """

	# Récupération des données
	n = data["n"] # n le nombre de clients	

	res_circuits = []
	for t =1:data["l"]
		circuits = Vector{Int64}[]
		for i in 2:n+1
			# On sait que tout chemin commence par le sommet 0 (centre de dépôt)
			if x[1,i,t] > 0
				passeDepot,circuit=find_tournee(n,x[:,:,t],i,1)
				push!(circuits, circuit)

			end
		end
		push!(res_circuits,circuits)
	end
	return res_circuits
end

function detect_sous_tour2()

end


#---------------------------------------------------------------
#PRP 2
#---------------------------------------------------------------


function PDI2_Exact_to_Circuit(data, x) 
    """
	Permet de représenter la solution renvoyé par le PL sous forme de liste de liste (représentant les circuits, dans l'ordre de passage)
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	x les valeurs des variables binaires xij qui vallent 1 si un véhicule utilise l’arc (i, j) 
	affichage qui vaut true si on veut afficher le graphe
    """

	# Récupération des données
	n = data["n"] # n le nombre de clients	

	res_circuits = [] #doit avoir data["l"] liste de circuits, dont chaque circuit contient <=k listes
	for t =1:data["l"]
		circuits = []
		for k =1:data["k"]
			#une voiture est utilisee
			if(sum(x[:,:,k,t])>1) #car x[:,1,k,t]+[1,:,k,t]==2
				#trouver le premier revendeur a visiter
				
				for i=2:n+1
					if(x[1,i,k,t]>0)
						passeDepot,circuit=find_tournee(n,x[:,:,k,t],i,1)
						push!(circuits,circuit)
						break
					end
				end
			end
		end
		push!(res_circuits,circuits)
	end
	return res_circuits
end



#---------------------------------------------------------------
#COMMUN SOUS TOURS, BigM
#---------------------------------------------------------------

function find_tournee(n,xkt,i,circuit)
	"""
	circuit: un chiffre, le point de depart du circuit
	n : nb de revendeur
	i : le point de depart suivant
	xkt : matrice x de la voiture k, a l instant t
	donne le circuit de la voirture k a l instant t
	"""
	depart=i
	tournee=[circuit]
	findDepot=false
	dejaPasse=Set() #stock les noeuds deja passes
	push!(dejaPasse,circuit)
	while(true)
		push!(tournee,depart)
		for j =1:n+1
			if(xkt[depart,j]>0)
				depart=j
				break
			end
		end
		if(depart==1)
			findDepot=true
			break
		elsif(depart in dejaPass)
			break
		end
	end
	return findDepot,tournee
end


function detect_sous_tour(n,xkt)
	"""
	n: nb de revendeurs
	xkt : matrice x[:,:,k,t]
	detecter un sous tour d une voiture k, a l instant t
	"""
	for i =1:n+1
		for j =1:n+1
			if(xkt[i,j]>1)
				return find_tournee(n,xkt,j,i)
			end
		end

	end

end


function Mt(data,t)
	s=0
	for T in [t : data["l"];]
		for i in [1 : data["n"];]
			s=s+data["d"][i,T]
		end
	end
	return min(data["C"],s)
end

function Mit_til(data,i,t)
	s=0
	for j in [t : data["l"];]
		s=s+data["d"][i,j]
	end
	return min(data["L"][i],data["Q"],s)
end



