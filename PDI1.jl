using JuMP
using CPLEX
using CPUTime
include("./tools.jl")
include("./toolsPDI.jl")


#--------------------------
#Principe d'indice :
# i \in N = {1,..,n,n+1}
# i \in Nc = {2,...,n+1}
#d : n+1,n+1
# L0 : n, i dans Nc
#
#
#
#--------------------------


function PLNE_PDI1(data,branch=false)
	"""
	data : dictionnaire de donnees
	branch : boolean, ajoute ou enleve les contraintes d elimination de sous tours
	"""
	lci=matrix_cout(data) #l'indice 1 est le centre de depot
						  # le cout du client i et j dans lci : ldi[i+1] et ldi[j+1]

	# Création d'un modèle. Ce modèle fera l'interface avec le solveur GLPK
	m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_EPINT" => 1e-15 ))
	#set_optimizer_attribute(m,"CPLEX_PARAM_TLIM",
	# Création des variables
	@variable(m, 0 <= p[1:data["l"]]) #quantite de production a la periode t
	@variable(m, y[1:data["l"]],Bin)
	@variable(m, z[1:data["n"], 1:data["l"]],Bin) # z_it=1 si le revendeur i est visite a la periode t, i in Nc
	@variable(m, 0 <= z0[1:data["l"]] <= data["k"],Int) #le nombre de vehicule qui part du depot noeud 0 a la periode t
	@variable(m, 0 <= w[1:data["n"] , 1: data["l"]] )
	
	@variable(m,0<=I[1:data["n"]+1,1:data["l"]])#inventaire du noeud i a la fin de la periode t
	@variable(m,0 <= q[1:data["n"]+1,1:data["l"]]) #quantite delivree au revendeur i a la periode tincl
	@variable(m, x[1:data["n"]+1,1:data["n"]+1,1:data["l"]],Bin) #x_ijt=1 si un vehicule voyage directement du noeud i au noeud j a l instant t
	
	@objective(m, Min, sum(data["u"]*p + data["f"]*y + sum( data["h"] .*I,dims=1 )[1,:] + sum(sum(lci .*x, dims=1), dims=2)[1,1,:] )) 
	
	#contraintes
	for t in [1:data["l"];]
		if(t==1)
			@constraint(m,data["L0"][1]+p[t] == sum(q[:,t])-q[1,t]+I[1,t] ) #(2) #t-1==0
			for i in [2: data["n"]+1;] 
				@constraint(m,data["L0"][i]+q[i,t] == data["d"][i-1,t] + I[i,t]) #(3) cas t-1==0
				@constraint(m,data["L0"][i]+q[i,t] <= data["L"][i]) #(6) #cas ou t-1==0
			end
		else 
			@constraint(m,I[1,t-1] + p[t] == sum(q[:, t])-q[1,t]+I[1,t]) #(2) cas t-1>0
		end
		@constraint(m,p[t]<=Mt(data,t)*y[t]) #(4)
		@constraint(m,I[1,t] <= data["L"][1]) #(5)
		@constraint(m,z0[t]<=data["k"]) #(10) 
		
		for i in [1: data["n"]+1;]
			if(i==1)
				@constraint(m,sum(x[i,k,t] for k =1:data["n"]+1 if i!=k)+ sum(x[k,i,t] for k =1:data["n"]+1 if k!=i)==2*z0[t]) #(9) cas ou i represente le centre de depot
				if(!branch)
					@constraint(m, w[i,t] <= data["Q"]*z0[t]) #(12)
				end
			elseif(i>1)
				if (t>1)
					@constraint(m,I[i,t-1]+q[i,t] == data["d"][i-1,t] + I[i,t]) #(3) t-1>0,i \in Nc
					@constraint(m,I[i,t-1]+q[i,t]<=data["L"][i]) #(6) #cas ou t-1>0, i dans Nc
				end
				@constraint(m,q[i,t]<=Mit_til(data,i-1,t)*z[i-1,t]) #(7) #i dans Nc, z[i-1,t] represente le revendeur i
				@constraint(m,sum(x[i,k,t] for k =1:data["n"]+1 if(k!=i))==z[i-1,t]) #(8)
				@constraint(m,sum(x[i,k,t] for k =1:data["n"]+1 if k!=i) + sum(x[k,i,t] for k =1:data["n"]+1 if k!=i)==2*z[i-1,t]) #(9)
				if(!branch)
					@constraint(m, w[i-1,t] <= data["Q"]*z[i-1,t]) #(12)
				end
			end
			if(!branch)
				for j in [2:data["n"]+1;]
					if(i!=j && i>1)
						@constraint(m,w[i-1,t]-w[j-1,t] >= q[i,t]-Mit_til(data,i-1,t)*(1-x[i,j,t])) #(11)
					end
				end
			end
		end
	end
	return m
end







