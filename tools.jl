using LinearAlgebra
using Graphs
using GraphPlot
using Cairo, Compose, Fontconfig
using Colors

function Read_file(filename)
	"""
	Lit le fichier d'instance et retourne un dictionnaire contenant toutes les données de l'instance en entrée
	Paramètres : 
    filename le chemin menant vers le fichier de l'instance à lire
	"""
	# Déclaration et initialisation
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
		#stock au depot non limité
		#capacite de production infinie
		d=Dict("type"=>typeA,"n" => n,"L0" => L0,"L" => L,"coord" => coord,"h" => h,"l" => l,"d" => d,"u" => u,"f" => f,"C" => C,"Q" => Q,"k" => k)
		return d
	elseif (typeA==2)
		d=Dict("type"=>typeA,"n" => n,"L0" => L0,"L" => L,"coord" => coord,"h" => h,"l" => l,"d" => d,"u" => u,"f" => f,"C" => C,"Q" => Q,"k" => k,"mc" => mc)
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

function WritePng_visualization_Graph(G, data, clients_t, qt, edge_colours, node_colours, coordx_t, coordy_t, filename)
	"""
	Permet d'enregistrer dans le répertoire courant un graphe passé en entrée en format Png
	Paramètres : 
	G un graphe (Graph)
    data le dictionnaire contenant les données de l'instance
	clients_t l'ensemble des noeuds à traiter pour ce pas de temps
	qt un tableau de demande pour chaque noeud, qt[i] représente la demande du noeud i pour ce pas de temps 
	edge_colours la couleur des arcs
	node_colours la couleur des noeuds
	coordx_t l'ensemble des coordonnées x des noeuds
	coordy_t l'ensemble des coordonnées y des noeuds
	filename le nom du fichier sous lequel on veut enregistrer l'image
	"""
	gp = gplot(G, coordx_t, coordy_t, nodelabel=[(clients_t[i], qt[i]) for i in 1:length(clients_t)], edgelinewidth = [200*length(clients_t) for e in edges(G)], edgestrokec=edge_colours, nodefillc=node_colours)
	draw(PNG(filename, 1.5*length(clients_t)cm, 1.5*length(clients_t)cm), gp) 
end

function WritePngGraph_Boites(data, q, t, circuits, filename)
	"""
	Permet, à partir d'un ensemble de boîtes, de créer le graphe associé et de l'enregistrer dans le répertoire courant
	Paramètres : 
    data le dictionnaire contenant les données de l'instance
	q un tableau [1:n, 1:l] de taille n*l, q[i,t] est la quantité produite pour le revendeur i à la période t
	t l'instant de temps à considérer
	circuits un ensemble de circuit représenté sous forme de liste de liste
	filename le nom du fichier sous lequel on veut enregistrer l'image
	"""
	# On créer l'ensemble de client à traiter (plus le centre de dépôt)
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

	# Position des noeuds définies par les données de l'instance
	coord = data["coord"]
	coordx_t = [0 for i in 1:length(clients_t)]
	coordy_t = [0 for i in 1:length(clients_t)]
	for i in clients_t
		coordx_t[dict[i]] = coord[i+1][1]
		coordy_t[dict[i]] = coord[i+1][2]
	end

	node_colours = [colorant"blue" for i in 1:length(clients_t)]
	dictEdge = Dict() # Dictionnaire qui associe à chaque arc sa couleur

	nb_circuits = length(circuits)
	num_circuit = 1

	couleurs = distinguishable_colors(nb_circuits, [RGB(1,1,1), RGB(0,0,0)], dropseed=true) # Pour avoir un ensemble de couleurs distinguables
	for circuit in circuits # Un des circuits au temps t
		
		node_colours[1] = RGBA(0,0,1) # Puisque dict[0]=1, on met donc le noeud du dépôt en bleu
		for (i, j) in zip(circuit[1:end-1], circuit[2:end])
			add_edge!(G, dict[i], dict[j]) 
			dictEdge[(dict[i], dict[j])] = couleurs[num_circuit] # On ajoute au dictionnaire de couleur la couleur de l'arc (i,j)
			node_colours[dict[j]] = couleurs[num_circuit] # La couleur du sommet j
		end

		add_edge!(G, dict[circuit[end]], 1) 
		dictEdge[(dict[circuit[end]], 1)] = couleurs[num_circuit]
		num_circuit += 1
	end

	# On défini les couleurs des arcs (dans l'ordre) grâce au dictionnaire précédemment établi
	edge_colours = []
	for (e_idx, e) in enumerate(edges(G))
        i = src(e)
        j = dst(e)
		push!(edge_colours, dictEdge[(i,j)])
	end

	# Création de qt, les q[i,t] pour le pas de temps t considéré et pour les clients considéré seulement 
	qt = [0.0] # Le dépôt n'a pas de demande
	for i in clients_t[2:length(clients_t)]
		push!(qt, q[i,t])
	end

	WritePng_visualization_Graph(G, data, clients_t, qt, edge_colours, node_colours, coordx_t, coordy_t, filename * "_" * string(t) * ".png")
end

# ----- Tests -----

# pathFileData = "PRP_instances/B_200_instance30.prp"

# data = Read_file(pathFileData)
# println(data)
# println(data["n"])

# println("matrice cout : ")
# println(matrix_cout(data))

# WritePngGraph_Boites(data, q, 2, boites, "graphe") # boites provient de VRP_Heuristiques.jl et q de LSP.jl