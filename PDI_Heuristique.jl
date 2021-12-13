using JuMP
using CPLEX
using DelimitedFiles

# include("tools.jl")
# include("LSP.jl")
# include("VRP_Heuristiques.jl")
# include("VRP.jl")

function PDI_heuristique(data, heuristique, filename)
    """
    Paramètres : 
    data le dictionnaire contenant les données de l'instance
	heuristique une chaine de caractère déterminant quelle heuristique on va utiliser (ou "PL" si on veut résoudre exactement avec le PL)
	filename le nom sous lequel on veut enregistrer les graphes (0 si on ne veut pas les enregistrer) 
    """

    nbIteMax = 3 #PB a adapter

    # Récupération des données 
    n = data["n"] # n le nombre de revendeurs (en comptant le fournisseur)
    l = data["l"] # l horizon de planification
    f = data["f"] # f coût fixe par période de production
    u = data["u"] # u coût unitaire
    h = data["h"] # h tableau de dimension n, h[i] exprime le coût de stockage unitaire du revendeur i
    c = matrix_cout(data) # c tableau de dimension n*n, c[i][j] exprime le coût de transport du nœud i au nœud j
    # Pour rappel, tous les indices de c et h doivent être augmentés de 1 par rapport aux indices de l'énoncé (en julia, les indices commencent à 1)

    # Initialisation des paramètres SC : 
    SC = Array{Float64}(undef,n,l)
    for i=1:n
        for t=1:l
            SC[i,t] = c[1,i+1] + c[i+1,1]
        end
    end

    bestSol = []
    CoutbestSol = 0
    nbIte = 0
    while nbIte < nbIteMax
        # Résolution du problème LSP avec prise en compte des coûts de visite
        p, y, I, q = PL_LSP(data, SC, false)
        # println("p=", p) 
        # println("y=", y)
        # println("I=", I)
        # println("q=", q)

        # Coûts engendrés par le LSP (stockage + production)
        CoutSolcourante = sum(u*p[t] + f*y[t] + sum(h[i+1]*I[i,t] for i = 1:n) for t = 1:l )
        #println("PB CoutSolcourante du LSP = ", CoutSolcourante)

        tournees = []
        # Résolution, pour chaque pas de temps, du VRP
        for t=1:l
            if sum(q[i,t] for i in 1:n) > 0 # Si une production est lancée pour des clients à la période t, sinon ne rien faire à ce temps
                # Résoudre VRP avec les paramètres : p[t] quantité produite à la période t, 
                # y[t] vaut 1 si une production est lancée à la période t, et 0 sinon,
                # I[:,t] la quantité en stock à la fin de la période t pour chaque revendeur i, et
                # q[:,t] la quantité produite pour chaque revendeur i à la période t
                # Avec VRP, on livre à chaque pas de temps t, q[i,t] unités pour chaque client i (qui prend en compte le stock, etc. dans LSP)
                
                # PB ce n'est pas pcq y[t] = 1 qu'on va forcément livrer aux clients (on peut décider de tout stocker)
                circuits_t = VRP_iteratif(data, q, t, heuristique)
                if circuits_t == 0
                    println("!! Problème dans la résolution du VRP !!")
                    return 0
                end

                push!(tournees, circuits_t) # On ajoute l'ensemble des tournées effectuées au temps t
                CoutSolcourante += Cout_Circuits(data, circuits_t)
                #println("PB CoutSolcourant après VRP = ", CoutSolcourante)
            else
                push!(tournees,[]) # Pour signifier qu'aucune tournée n'est effectuée au temps t
            end
        end

        if nbIte == 0 # Initialisation de la meilleure solution
            bestSol = [p, y, I, q, tournees] # Mémorisation de la meilleure solution
            CoutbestSol = CoutSolcourante
        elseif CoutbestSol > CoutSolcourante
            bestSol = [p, y, I, q, tournees] # Mémorisation de la meilleure solution
            CoutbestSol = CoutSolcourante
        end
        
        # Mise à jour des coûts SC
        # On a length(tournees) qui vaut l
        for t_tournee in 1:length(tournees) # t_tournee est le pas de temps considéré pour l'ensemble de tournée en cours d'étude
            elems_t = [] # Ensemble des revendeurs livrés en t
            for circuit in tournees[t_tournee] 
                for ind_elem in 1:length(circuit) 
                    push!(elems_t, circuit[ind_elem])
                    if ind_elem == 1 || circuit[ind_elem] == 0 # Seule la première partie du OU (||) suffit 
                        # Dans ce cas circuit[ind_elem] vaut 0 (le dépôt) et on ne veut pas de SC[i,t] pour i = 0, on passe
                        continue
                    elseif ind_elem == length(circuit)
                        # Le premier indice de SC va de 1 à n (sans le 0) donc pas +1 mais plutôt ne pas traiter qd circuit[ind_elem]=0
                        SC[circuit[ind_elem], t_tournee] = c[circuit[ind_elem-1]+1, circuit[ind_elem]+1] + c[circuit[ind_elem]+1, circuit[1]+1] - c[circuit[ind_elem-1]+1, circuit[1]+1]
               
                    else # Si ind_elem != 1 && ind_elem != length(circuit)
                        # println("PB c[circuit[length(circuit)]+1, circuit[ind_elem]+1]",c[circuit[length(circuit)]+1, circuit[ind_elem]+1])
                        # println("PB c[circuit[ind_elem]+1, circuit[ind_elem+1]+1] ",c[circuit[ind_elem]+1, circuit[ind_elem+1]+1])
                        # println("PB circuit[length(circuit)+1]", circuit[length(circuit)]+1)
                        # println("PB circuit[ind_elem+1]+1",circuit[ind_elem+1]+1)
                        # println("PB c[circuit[length(circuit)]+1, circuit[ind_elem+1]+1]",c[circuit[length(circuit)]+1, circuit[ind_elem+1]+1])
                        # println("PB circuit[",ind_elem,"]", circuit[ind_elem])
                        # println("PB SC[circuit[ind_elem], t_tournee]", SC[circuit[ind_elem], t_tournee])
                        SC[circuit[ind_elem], t_tournee] = c[circuit[ind_elem-1]+1, circuit[ind_elem]+1] + c[circuit[ind_elem]+1, circuit[ind_elem+1]+1] - c[circuit[ind_elem-1]+1, circuit[ind_elem+1]+1]
                    end 
                end
            end

            # On parcourt tous les noeuds et on màj les SCit des noeuds i ne faisant pas partie de elems_t, càd n'étant pas livré au temps t            
            for i in 1:n
                if !(i in elems_t)
                    # Si tournees[t_tournee] = [] c'est qu'aucun revendeur n'a été livré au temps t
                    if tournees[t_tournee] == []
                        # La valeur du coût de l’insertion la plus économique dans ce cas est simplement le coût du circuit passant par 0 et par i                
                        SC[i, t_tournee]= c[1,i+1] + c[i+1,1]
                    else
                        circuit = tournees[t_tournee][1]

                        # Initialisation du coût
                        if length(circuit) == 1 # C'est à dire que le circuit n'est composé que du dépôt 
                            # Cout_Min sera le coûts d'insertion de i le plus économique
                            Cout_Min = c[circuit[1]+1, i+1] + c[i+1, circuit[1]+1] # - c[circuit[1]+1, circuit[1]+1] vaut 0 ici 
                        else
                            Cout_Min = c[circuit[1]+1, i+1] + c[i+1, circuit[2]+1] - c[circuit[1]+1, circuit[2]+1] # Cout_Min sera le coûts d'insertion de i le plus économique
                        end
                        
                        for circuit in tournees[t_tournee] 
                            for ind_elem in 1:length(circuit)
                                if ind_elem == length(circuit)
                                    Cout_Courrant = c[circuit[ind_elem]+1, i+1] + c[i+1, circuit[1]+1] - c[circuit[ind_elem]+1, circuit[1]+1]
                        
                                elseif ind_elem != length(circuit)
                                    Cout_Courrant = c[circuit[ind_elem]+1, i+1] + c[i+1, circuit[ind_elem+1]+1] - c[circuit[ind_elem]+1, circuit[ind_elem+1]+1]
                                end

                                if Cout_Courrant < Cout_Min
                                    Cout_Min = Cout_Courrant 
                                end
                            end
                        end
                        
                        SC[i, t_tournee]= Cout_Min
                    end
                end
            end
        end

        nbIte += 1
    end

    if filename != 0
        # Une solution est de la forme [p, y, I, q, tournees]
        p, y, I, q, tournees = bestSol

        # pt (Variable de production): quantité produite à la période t.
        # yt (Variable de lancement): variable binaire qui vaut 1 si une production est lancée à la période t, et 0 sinon. (se déduit de pt)
        # Iit (Variable de stockage): quantité en stock à la fin de la période t pour i. (se déduit du stock de départ, de la demande en i et de qit)
        # qit (Variable d’approvisionnement): quantité produite pour le revendeur i à la période t.
        # tournees ensemble des tournées (càd ensemble de tournées pour tous les pas de temps).

        # On représente donc une solution par un graphe décrivant l'ensemble des tournées pour chaque pas de temps et 
        # avec la quantité produite pour le revendeur i à la période t (q[i,t]) indiquée sur chaque noeud i
        for t in 1:length(tournees)
            if tournees[t] != []
                WritePngGraph_Boites(data, q, t, tournees[t], filename)
            end
        end
        println("Quantité à produire à chaque période : ", p)
    end

    return bestSol, CoutbestSol
end

# ----- Tests -----

#pathFileData = "PRP_instances/B_200_instance30.prp"
#pathFileData = "PRP_instances/A_014_ABS1_15_1.prp"
# pathFileData = "PRP_instances/B_100_instance30.prp" #693092 pour nbite=2 #676261 pour nbite = 3,4,5 et 15
# data = Read_file(pathFileData)

# bestSol, CoutbestSol = PDI_heuristique(data, "CW", 0)
# println("Solution de coût ", CoutbestSol)
# PB il faut aussi donner les pt pou caractériser une solution !


# Récupération des coût et des temps d'éxecution
# Pour TSP approché (plus proche voisin parmi les restants)/TSP exact #PB
# Pour VRP PLNE, Bin Packing, Clark Wright, Sectorielle ("PL", "BP", "CW", "Sec")
# Dans VRP_iteratif pour le changement de boîte d'un client, soit Cout_heur_Boite soit MinCout_Ajout
# Si length(boites) < k, ajout de boîtes vides ou couper les boîtes en 2
heurs = ["BP", "CW", "Sec"] #, "PL"] PB PL trop long je pense

function Calculs_Temps_Cout(files, heurs)
    resTempsCout = Dict()
    for file in files
        for heur in heurs        
            pathFileData = "PRP_instances/" * file * ".prp"
            data = Read_file(pathFileData)

            time1 = Sys.time()
            bestSol, CoutbestSol = PDI_heuristique(data, heur, 0)
            time2 = Sys.time()

            Tdiff = time2 - time1
            resTempsCout[(heur, file)] = (Tdiff, CoutbestSol)
        end
    end
    return resTempsCout
end

# files1 = ["A3", "A4", "A5", "A6", "A10", "A62"]
# resTempsCout1 = Calculs_Temps_Cout(files1, heurs)
# println(resTempsCout1)
# # Enregistrer les données dans un fichier texte
# writedlm("Temps_Cout_1.txt", resTempsCout1)

files2 = ["A_014_ABS96_15_5", "A_050_ABS19_50_5", "A_014_ABS11_15_1", "A_014_ABS21_15_1", "A_014_ABS31_15_1", "A_014_ABS41_15_1", "A_014_ABS51_15_1", "A_014_ABS61_15_1", "A_014_ABS71_15_1", "A_014_ABS81_15_1", "A_014_ABS91_15_1", "A_100_ABS92_100_1"]
resTempsCout2 = Calculs_Temps_Cout(files2, heurs)
println(resTempsCout2)
# Enregistrer les données dans un fichier texte
writedlm("Temps_Cout_2.txt", resTempsCout2)

files3 = ["A_050_ABS13_50_3","B_050_instance2","B_050_instance7","B_050_instance14","B_050_instance30","B_050_instance22", "B_100_instance20", "A_100_ABS1_100_1"]
resTempsCout3 = Calculs_Temps_Cout(files3, heurs)
println(resTempsCout3)
# Enregistrer les données dans un fichier texte
writedlm("Temps_Cout_3.txt", resTempsCout3)

files4 = ["B_100_instance1","B_100_instance6","B_100_instance16","B_100_instance30","A_014_ABS71_15_4","A_014_ABS80_15_4","A_014_ABS94_15_5"]
resTempsCout4 = Calculs_Temps_Cout(files4, heurs)
println(resTempsCout4)
# Enregistrer les données dans un fichier texte
writedlm("Temps_Cout_4.txt", resTempsCout4)

files5 = ["A_050_ABS1_50_2", "A_050_ABS77_50_3", "A_050_ABS78_50_2", "A_100_ABS12_100_4", "A_100_ABS20_100_3", "A_100_ABS28_100_3"]
resTempsCout5 = Calculs_Temps_Cout(files5, heurs)
println(resTempsCout5)
# Enregistrer les données dans un fichier texte
writedlm("Temps_Cout_5.txt", resTempsCout5)