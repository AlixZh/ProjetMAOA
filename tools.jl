using LinearAlgebra
using Graphs
using GraphPlot
#using SimpleWeightedGraphs # PB SimpleWeightedDiGraphs ?
using Cairo, Compose, Fontconfig
using Colors

function Read_file(filename)
	"""
	Lit le fichier d'instance et retourne un dictionnaire contenant toutes les données de l'instance en entrée
	Paramètres : 
    filename le chemin menant vers le fichier de l'instance à lire
	"""
	#declaration et initialisation
	ligne_demande=false 
	typeA=0
	n=0 # nb clients
	l=0 # nb horizon
	u=0 # cout unitaire de production
	f=0 # cout fixe de production
	C=1e+10 #capacite de production
	Q=0 #capacite de maximale de chaque vehicule
	k=0 # nb de vehicule   
	L=[] #capacité de stockage maximale
	L0=[] #capacite de stockage maximale au centre de depot
	coord=[] #coordonnees des positions des clients et du centre de depot
	d=[] #demande du client i au temps t
	h=[] # coût de stockage unitaire
	mc=0
    open(filename) do fichier
		for (i,line) in enumerate(eachline(fichier))
			x = split(line," ") 
			if(x[1]=="Type")
				typeA=Int(parse(Float64,x[2]))
			elseif(x[1]=="n")       # A line beginning with a 'n' gives the graph size -1
				n = parse(Int,x[2]) 
				#g = SimpleGraph(n+1)  # creation of a undirected graph with n+1 nodes 
				L0 = zeros(n+1)
				L = zeros(n+1)
				coord = Vector{Tuple}(undef,n+1)
				h = zeros(n+1) # coût de stockage unitaire
			elseif(x[1] == "l") # A line beginning with a 'l' gives the edges
				l = Int(parse(Float64,x[2])) # nb horizon
				d=zeros(n,l)
			elseif(x[1]=="u")
				u = parse(Int,x[2]) # cout unitaire de production
			elseif(x[1]=="f")
				f = parse(Int,x[2]) # cout fixe de production
			elseif(x[1]=="C")
				C = Int(parse(Float64,x[2])) #capacite de production
			elseif(x[1]=="Q")
				Q = parse(Int,x[2]) # capacite maximale de chaque vehicule
			elseif(x[1]=="k")
				k = parse(Int,x[2]) # nb de vehicule               
			elseif(x[1]=="d")
				ligne_demande=true
			elseif(x[1]=="mc")
				mc=parse(Int,x[2])			
			elseif(ligne_demande)
				x1=parse(Int,x[1])
				for t in [1:l;]
					d[x1,t] = parse(Float64,x[1+t])
				end
			else
				x1 = parse(Int,x[1])
				L[x1+1] = parse(Float64,x[8])
				L0[x1+1] = parse(Float64,x[10])
				h[x1+1] = parse(Float64,x[6])
				coord[x1+1] = (parse(Float64,x[2]) ,parse(Float64,x[3]))   
			end
		end
	end
	if(typeA==1)
		#stock au depot non limite
		#capacite de production infinie
		d=Dict("type"=>typeA,"n" => n,"L0" => L0,"L" => L,"coord" => coord,"h" => h,"l" => l,"d" => d,"u" => u,"f" => f,"C" => C,"Q" => Q,"k" => k)
		#return n,L0,L,coord,h,l,d,u,f,C,Q,k
		return d
	elseif (typeA==2)
		d=Dict("type"=>typeA,"n" => n,"L0" => L0,"L" => L,"coord" => coord,"h" => h,"l" => l,"d" => d,"u" => u,"f" => f,"C" => C,"Q" => Q,"k" => k,"mc" => mc)
		#return n,L0,L,coord,h,l,d,u,f,C,Q,k,mc
		return d
	else
		println("probleme fichier instance : ",filename)
		return

	end
end

function coutA(xi,xj)
	"""
	Retourne le coût (pour les instances de type A) entre les deux points en entrée
	Paramètres : 
	xi un tuple de deux valeurs contenant les coordonnées du client i (xi[1], xi[2])
	xj un tuple de deux valeurs contenant les coordonnées du client j (xj[1], xj[2])
	"""
	return floor(sqrt( (xi[1]-xj[1])*(xi[1]-xj[1]) + (xi[2]-xj[2])*(xi[2]-xj[2]) ) + (1/2) )
end 


function coutB(xi,xj,mc)
	"""
	Retourne le coût (pour les instances de type B) entre les deux points en entrée
	Paramètres : 
	xi un tuple de deux valeurs contenant les coordonnées du client i (xi[1], xi[2])
	xj un tuple de deux valeurs contenant les coordonnées du client j (xj[1], xj[2])
	mc une constante multiplicative (donnée dans l'instance à considérer)
	"""
	return mc*sqrt( (xi[1]-xj[1])*(xi[1]-xj[1]) + (xi[2]-xj[2])*(xi[2]-xj[2]) )
end 

function matrix_cout(data)
	"""
	Retourne la matrice de coût associé aux données de l'instance en entrée
	Paramètres : 
    data le dictionnaire contenant les données de l'instance
	"""
	c=zeros(data["n"]+1,data["n"]+1)
	for i in [1:data["n"]+1;] 
		for j in [1:data["n"]+1;] 
			if(data["type"] == 1)
				c[i,j] = coutA(data["coord"][i], data["coord"][j])
			elseif(data["type"] == 2)
				c[i,j] = coutB(data["coord"][i], data["coord"][j], data["mc"])
			else 
				println("Probleme de type dans le calcul des couts")
				return c
			end
		end
	end
	return c
end

function WritePdf_visualization_Graph(G, filename)
	"""
	Permet d'enregistrer dans le répertoire courant un graphe passé en entrée en format Pdf
	Paramètres : 
	G un graphe (Graph)
	filename le nom du fichier sous lequel on veut enregistrer l'image
	"""
	filename_splitted_in_two_parts = split(filename,".") # split to remove the file extension
	filename_with_pdf_as_extension= filename_splitted_in_two_parts[1]*".pdf"
	
	# save to pdf
	draw(PDF(filename_with_pdf_as_extension, 16cm, 16cm), gplot(G, nodelabel = 0:nv(G)))
end

# PB
# function WritePng_visualization_Graph(G, data, qt, edge_colours, filename)
# 	"""
# 	Permet d'enregistrer dans le répertoire courant un graphe passé en entrée en format Png
# 	Paramètres : 
# 	G un graphe (Graph)
#     data le dictionnaire contenant les données de l'instance
# 	edge_colours la couleur des arcs
# 	filename le nom du fichier sous lequel on veut enregistrer l'image
# 	"""
# 	draw(PNG(filename, 16cm, 16cm), gplot(G, nodelabel=[(i, qt[i+1]) for i in 0:data["n"]]))#, edgestrokec=edge_colours)) #PB ou ajouter au nodelabel q[i,t]
# end

# function WritePngGraph_Boites(data, q, t, circuits, filename) #PB affichage a partir de PDI_heuristique
# 	"""
# 	Permet, à partir d'un ensemble de boîtes, de créer le graphe associé et de l'enregistrer dans le répertoire courant
# 	Paramètres : 
#     data le dictionnaire contenant les données de l'instance
# 	q un tableau [1:n, 1:l] de taille n*l, q[i,t] est la quantité produite pour le revendeur i à la période t
# 	t l'instant de temps à considérer
# 	circuits un ensemble de circuit représenté sous forme de liste de liste
# 	filename le nom du fichier sous lequel on veut enregistrer l'image
# 	"""

# 	#G = SimpleWeightedGraph(data["n"] + 1) #PB avec des couts sur les arcs (q[i,t] transporté)

# 	G = Graph(data["n"] + 1)

# 	edge_colours = []

# 	nb_circuits = length(circuits)
# 	num_circuit = 0
# 	for circuit in circuits # Un des circuit au temps t
# 		# sommeqit = 0
# 		# for i in circuit
# 		# 	if i != 0 # Le dépôt n'a pas de demande
# 		# 		sommeqit += q[i,t] # L'ensemble de charge à transporter par le camion sur ce circuit (inférieur à Q)
# 		# 	end
# 		# end

# 		couleur = RGBA(0, num_circuit/nb_circuits, num_circuit/nb_circuits) # La couleur pour ce circuit

# 		for (i, j) in zip(circuit[1:end-1], circuit[2:end])
# 			add_edge!(G, i+1, j+1) #, sommeqit)
# 			push!(edge_colours, couleur)

# 			if j == 0 #PB
# 				println("PB j=", j)
# 				println("PB circuit = ", circuit)
# 			end
			
# 			# if j != 0 # Le dépôt n'a pas de demande # PB normalement pas besoin puisqu'aucun arc allant vers 0 sauf le dernier arc qui n'est pas compté ici
# 			# 	sommeqit -= q[j,t] # On dépose q[j,t] au sommet j
# 			# end
# 		end

# 		add_edge!(G, circuit[end]+1, circuit[1]+1) #, sommeqit) # Normalement sommeqit vaut 0 ici puisque c'est l'arc qui revient au dépôt
# 		push!(edge_colours, couleur)
# 		num_circuit += 1
# 	end

# 	# Création de qt
# 	qt = [0]

# 	for i in 1:data["n"]
# 		push!(qt, q[i,t])
# 	end

# 	WritePng_visualization_Graph(G, data, qt, edge_colours, filename * "_" * string(t))
# end

function WritePng_visualization_Graph(G, data, clients_t, qt, edge_colours, filename)
	"""
	Permet d'enregistrer dans le répertoire courant un graphe passé en entrée en format Png
	Paramètres : 
	G un graphe (Graph)
    data le dictionnaire contenant les données de l'instance
	clients_t l'ensemble des noeuds à traiter pour ce pas de temps
	qt un tableau de demande pour chaque noeud, qt[i] représente la demande du noeud i pour ce pas de temps 
	edge_colours la couleur des arcs
	filename le nom du fichier sous lequel on veut enregistrer l'image
	"""
	gp = gplot(G, layout = stressmajorize_layout, nodelabel=[(clients_t[i], qt[i]) for i in 1:length(clients_t)], edgelinewidth = [10000*length(clients_t) for e in edges(G)], edgestrokec=edge_colours, nodefillc=edge_colours)
	draw(PNG(filename, 1.5*length(clients_t)cm, 1.5*length(clients_t)cm), gp) 
end

function WritePngGraph_Boites(data, q, t, circuits, filename) #PB affichage a partir de PDI_heuristique
	"""
	Permet, à partir d'un ensemble de boîtes, de créer le graphe associé et de l'enregistrer dans le répertoire courant
	Paramètres : 
    data le dictionnaire contenant les données de l'instance
	q un tableau [1:n, 1:l] de taille n*l, q[i,t] est la quantité produite pour le revendeur i à la période t
	t l'instant de temps à considérer
	circuits un ensemble de circuit représenté sous forme de liste de liste
	filename le nom du fichier sous lequel on veut enregistrer l'image
	"""

	#G = SimpleWeightedGraph(data["n"] + 1) #PB avec des couts sur les arcs (q[i,t] transporté)

	# On créer l'ensemble de client à traiter
	clients_t = [0]
	for i in 1:data["n"]
		if q[i,t] != 0 # On ne traite que les clients qui ont une demande pour ce pas de temps
			push!(clients_t, i)
		end
	end

	dict = Dict() # Dictionnaire qui associe à chaque élément son indice dans clients_t
	for ind in 1:length(clients_t)
		if !(clients_t[ind] in keys(dict))
			dict[clients_t[ind]] = ind
		end
	end

	G = DiGraph(length(clients_t))

	# Position des noeuds #PB pos des noeuds definie
	coord = data["coord"]
	coords_t = []
	for i in clients_t
		push!(coords_t, coord[i+1])
	end

	edge_colours = []

	nb_circuits = length(circuits)
	println("PB NB CIRCUIT = ", nb_circuits)
	num_circuit = 1
	# r = 1
	# g = 0
	# b = 0

	couleurs = distinguishable_colors(nb_circuits, colorant"blue")
	for circuit in circuits # Un des circuit au temps t
		# sommeqit = 0
		# for i in circuit
		# 	if i != 0 # Le dépôt n'a pas de demande
		# 		sommeqit += q[i,t] # L'ensemble de charge à transporter par le camion sur ce circuit (inférieur à Q)
		# 	end
		# end

		# couleur = RGBA(r, g, b) # La couleur pour ce circuit

		# if num_circuit == 0
		# 	r = 0
		# 	b = 1
		# 	g = 0
		# end
		# if num_circuit == 1
		# 	r = 0
		# 	b = 0
		# 	g = 1
		# end
		# if num_circuit%3 == 0
		# 	r = r + 0.5
		# 	if r>1
		# 		r = 0
		# 	end
		# elseif num_circuit%3 == 1
		# 	b = b + 0.5
		# 	if b>1
		# 		b = 0
		# 	end
		# elseif num_circuit%3 == 2
		# 	g = g + 0.5
		# 	if g>1
		# 		g = 0
		# 	end
		# end

		for (i, j) in zip(circuit[1:end-1], circuit[2:end])
			add_edge!(G, dict[i], dict[j]) #, sommeqit)
			#println("PB Ajout arc entre ", dict[i], " càd ", i, " et ", dict[j], " càd ", j)
			push!(edge_colours, couleurs[num_circuit])

			if j == 0 #PB
				println("PB j=", j)
				println("PB circuit = ", circuit)
			end
			
			# if j != 0 # Le dépôt n'a pas de demande # PB normalement pas besoin puisqu'aucun arc allant vers 0 sauf le dernier arc qui n'est pas compté ici
			# 	sommeqit -= q[j,t] # On dépose q[j,t] au sommet j
			# end
		end

		add_edge!(G, dict[circuit[end]], 1) #, sommeqit) # Normalement sommeqit vaut 0 ici puisque c'est l'arc qui revient au dépôt
		#println("PB Ajout arc entre ", dict[circuit[end]], " càd ", circuit[end], " et 1 ")
		push!(edge_colours, couleurs[num_circuit])
		num_circuit += 1
	end

	# Création de qt
	qt = [0]

	for i in clients_t[2:length(clients_t)]
		push!(qt, q[i,t])
	end

	WritePng_visualization_Graph(G, data, clients_t, qt, edge_colours, filename * "_" * string(t))
end

#Tests
pathFileData = "PRP_instances/B_200_instance30.prp"
# pathFileData = "PRP_instances/A_014_ABS2_15_1.prp"

data = Read_file(pathFileData)
# println(data)
# println(data["n"])

# println("matrice cout : ")
# println(matrix_cout(data))

WritePngGraph_Boites(data, q, 2, boites, "graphe") # boites provient de VRP_Heuristiques.jl et q de LSP.jl












