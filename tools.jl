using LinearAlgebra
using Graphs
using GraphPlot
using Cairo, Compose, Fontconfig

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
	Permet d'enregistrer un graphe passé en entrée en format Pdf
	Paramètres : 
	G un graphe (Graph)
	filename le nom du fichier sous lequel on veut enregistrer l'image
	"""
	filename_splitted_in_two_parts = split(filename,".") # split to remove the file extension
	filename_with_pdf_as_extension= filename_splitted_in_two_parts[1]*".pdf"
	
	# save to pdf
	draw(PDF(filename_with_pdf_as_extension, 16cm, 16cm), gplot(G, nodelabel = 0:nv(G)))
end

function WritePng_visualization_Graph(G, data, filename)
	"""
	Permet d'enregistrer un graphe passé en entrée en format Png
	Paramètres : 
	G un graphe (Graph)
    data le dictionnaire contenant les données de l'instance
	filename le nom du fichier sous lequel on veut enregistrer l'image
	"""
	draw(PNG(filename, 16cm, 16cm), gplot(G, nodelabel=0:data["n"]))
end

function WritePngGraph_Boites(data, circuits, filename)
	"""
	Permet, à partir d'un ensemble de boîtes, de créer le graphe associé et de l'enregistrer
	Paramètres : 
    data le dictionnaire contenant les données de l'instance
	circuits un ensemble de circuit représenté sous forme de liste de liste
	filename le nom du fichier sous lequel on veut enregistrer l'image
	"""

	G = Graph(data["n"] + 1)

	for circuit in circuits
		for (i, j) in zip(circuit[1:end-1], circuit[2:end])
			add_edge!(g, i+1, j+1)
		end

		add_edge!(G, circuit[end]+1, circuit[1]+1)
	end

	WritePng_visualization_Graph(G, data, filename)
end

#Tests
pathFileData = "PRP_instances/B_200_instance30.prp"
# pathFileData = "PRP_instances/A_014_ABS1_15_1.prp"

data = Read_file(pathFileData)
# println(data)
# println(data["n"])

# println("matrice cout : ")
# println(matrix_cout(data))

WritePngGraph_Boites(data, boites, "graphe") #boites provient de VRP_Heuristiques













