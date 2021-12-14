using Plots
using DelimitedFiles

include("BC.jl")

function Dessine_courbe(tabX,tabY)
	nom="Temps d'execution en fonction de la taille de l'instance"
    println("Création du fichier de courbes: ", nom)

    Plots.plot(tabX,tabY,label="instanceA",title=nom)
    xlabel!(nom)
    Plots.savefig(nom*".pdf")  # Ecrit la courbe créée à la ligne précédente dans un fichier .pdf
	
end   


function execution(numAlgo,sortiename,path="./instancesTests/",pathDest="./Sortie/")
	"""
	tracer les courbes et afficher le temps d execution et la valeur opt
	path : chemin des fichiers a tester
	"""
	res=Dict()
	listName=[]
	
  
    
	for namefile in readdir(path)
		data = Read_file(path*namefile)
		if(numAlgo==1)
			trouve,opt,tps,x,y,q,p,gap = BC1(data)
		elseif(numAlgo==2)
			trouve,opt,tps,x,y,q,p,gap = BC2(data)
		end
		
		if(trouve)
			println()
			println("Nom du fichier : ", namefile)
			res[namefile]=Dict()
			push!(listName,namefile)
			res[namefile]["name"]=namefile
			res[namefile]["opt"]=opt
			res[namefile]["temps"]=tps
			res[namefile]["x"]=x
			res[namefile]["q"]=q
			res[namefile]["y"]=y
			res[namefile]["data"]=data
			res[namefile]["p"]=p
			res[namefile]["gap"]=gap
			if(numAlgo==1)
				res[namefile]["ens_circuit"]=PDI1_Exact_to_Circuit(data, x) 
			elseif(numAlgo==2)
				res[namefile]["ens_circuit"]=PDI2_Exact_to_Circuit(data, x) 
			end
			myq=zeros(data["n"],data["l"])
			for i=2:data["n"]+1
				for j=1:data["l"]
					myq[i-1,j]=q[i,j]
				end
			
			end
			for t=1:data["l"]
				if(length(res[namefile]["ens_circuit"][t])>=1)
					WritePngGraph_Boites(data,myq , t, res[namefile]["ens_circuit"][t], pathDest*"Boite_"*namefile)
				end
			end
			
			println("Quantité à produire à chaque période : ", p)
			println()
		end
	end
	open(pathDest*sortiename, "w") do io
    	writedlm(io, res)
    	writedlm(io,listName)
    end
	return listName,res
end
