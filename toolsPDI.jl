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
				for i =1 : length(circuit)
					circuit[i]=circuit[i]-1
				end
				push!(circuits, circuit)

			end
		end
		push!(res_circuits,circuits)
	end
	return res_circuits
end



function detect_sous_tour1(n,xkt)
	"""
	n: nb de revendeurs
	xkt : matrice x[:,:,k,t], contient peut etre plusieurs tournees
	detecter un sous tour d une voiture k, a l instant t
	"""
	sousTours=Set()
	dejaVisite=Set()
	trouveDepartDepot=true # =false si trouver au moins un sous tour
	for i =2:n+1
		if(!(i in dejaVisite) ) #si deja visite, appartient deja a un sous tours
			for j =2:n+1
				if(xkt[i,j]>=1)
					passeDepot,tournee=find_tournee(n,xkt,j,i)
					if(!passeDepot && !(1 in tournee)) #detection d un sous tours
						trouveDepartDepot=false
						push!(sousTours, Set(tournee))
						union!(dejaVisite,tournee)

					end
				end
			end
		end

	end
	#print(sousTours)
	return trouveDepartDepot,sousTours
end


function get_param_opt1(data,m)
	"""
	m : model
	data : donnees de l instance
	"""
	xopt=zeros(Int,data["n"]+1,data["n"]+1,data["l"])
	yopt=zeros(Int,data["l"])
	qopt=zeros((data["n"]+1,data["l"]))
	popt=zeros(data["l"])
	for t =1:data["l"]
		yopt[t]=value(variable_by_name(m,"y[$t]"))
		popt[t] = value(variable_by_name(m,"p[$t]"))
		for i=1:data["n"]+1
			qopt[i,t]=value(variable_by_name(m,"q[$i,$t]"))
			for j=1:data["n"]+1
				xopt[i,j,t]=round(value(variable_by_name(m,"x[$i,$j,$t]")))
			end
		end
	end
	return xopt,yopt,qopt,popt
end



#---------------------------------------------------------------
#PRP 2
#---------------------------------------------------------------


function PDI2_Exact_to_Circuit(data, x) 
    """
	Permet de représenter la solution renvoyé par le PL sous forme de liste de liste (représentant les circuits, dans l'ordre de passage)
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	x les valeurs des variables binaires xij qui vaut 1 si un véhicule utilise l’arc (i, j) 
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




function detect_sous_tour2(n,xkt)
	"""
	n: nb de revendeurs
	xkt : matrice x[:,:,k,t]
	detecter un sous tour d une voiture k, a l instant t
	xkt contient une tournee de la voiture k
	"""
	for i =2:n+1
		for j =1:n+1
			if(xkt[i,j]>0)
				return find_tournee(n,xkt,j,i)
			end
		end

	end

end

function get_param_opt2(data,m)
	"""
	m : model
	data : donnees de l instance
	"""
	xopt=zeros(Int,data["n"]+1,data["n"]+1,data["k"],data["l"])
	yopt=zeros(Int,data["l"])
	qopt=zeros((data["n"]+1,data["k"],data["l"]))
	popt=zeros(data["l"])
	for t =1:data["l"]
		yopt[t]=value(variable_by_name(m,"y[$t]"))
		popt[t]=value(variable_by_name(m,"p[$t]"))
		for k =1:data["k"]
			for i=1:data["n"]+1
				qopt[i,k,t]=value(variable_by_name(m,"q[$i,$k,$t]"))
				for j=1:data["n"]+1
					xopt[i,j,k,t]=value(variable_by_name(m,"x[$i,$j,$k,$t]"))
				end
			end
		end
	end
	return xopt,yopt,qopt,popt
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
	#println(xkt)
	while(!findDepot && !(depart in dejaPasse))
		push!(tournee,depart)
		push!(dejaPasse,depart)
		for j =1:n+1
			if(xkt[depart,j]>0)
				depart=j
				break
			end
		end
		if(depart==1)
			if(!(1 in tournee))
				push!(tournee,depart)
			end
			findDepot=true
		end
	end
	#println(tournee)
	return findDepot,tournee
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








