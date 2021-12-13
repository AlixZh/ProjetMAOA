using JuMP
using CPLEX
using CPUTime
using IterTools
using Combinatorics

include("tools.jl")
include("./toolsPDI.jl")


#--------------------------
#Principe d'indice :
# i \in N = {1,..,n,n+1}
# i \in Nc = {2,...,n+1}
#d : n,n
# L0 : n+1, i dans N
#L :n+1
#lci: n+1,n+1, i \in N
#--------------------------


function PLNE_PDI2(data,branch=false)
	"""
	data : dictionnaire de donnees
	branch : boolean, ajoute ou enleve les contraintes d elimination de sous tours
	"""
	lci=matrix_cout(data) #l'indice 1 est le centre de depot
						  # le cout du client i et j dans lci : ldi[i+1] et ldi[j+1]

	# Création d'un modèle. Ce modèle fera l'interface avec le solveur GLPK
	m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_EPINT" => 1e-15 ))
	# Création des variables
	@variable(m,y[1:data["l"]],Bin) #(32)
	@variable(m,z[1:data["n"]+1, 1:data["k"],1:data["l"]],Bin)
	@variable(m,x[1:data["n"]+1, 1:data["n"]+1, 1:data["k"],1:data["l"]],Bin)
	@variable(m, 0 <= p[1:data["l"]] ) #(31)
	@variable(m, 0 <= I[1:data["n"]+1, 1:data["l"]])
	@variable(m, 0 <= q[1:data["n"]+1, 1:data["k"], 1:data["l"]])

	@objective(m, Min, sum(data["u"]*p + data["f"]*y + sum( data["h"] .*I,dims=1 )[1,:] + sum(sum(lci .*sum(x,dims=3), dims=1), dims=2)[1,1,1,:] )) 
	#@objective(m, Min, sum(data["u"]*p + data["f"]*y + sum( data["h"] .*I,dims=1 )[1,:] + sum( lci[i,j]*sum(x[i,j,:,t] for k =1:data["k"]) for i =1:data[n+1]+1, for j =1:data["n"]+1 if (i!=j)) )) 
	

	#contraintes

	for t = 1:data["l"]
		if(t>1)
			@constraint(m, I[1,t-1]+p[t] == sum(q[i+1,k,t] for i=1:data["n"] for k=1:data["k"]) + I[1,t]) #21
		else
			@constraint(m, data["L0"][1]+p[t] == sum(q[i+1,k,t] for i=1:data["n"] for k=1:data["k"]) + I[1,t]) #21
		end
		@constraint(m,p[t]<=Mt(data,t)*y[t]) #23
		@constraint(m, I[1,t]<=data["L"][1]) #24

		for i =2: data["n"]+1
			
			@constraint(m,sum(z[i,:,t])<=1) #(27)
			if( t>1)
				@constraint(m,I[i,t-1]+sum(q[i,:,t]) <= data["L"][i]) #(25)
				@constraint(m,I[i,t-1]+sum(q[i,:,t])==data["d"][i-1,t]+I[i,t]) #22
			else
				@constraint(m,data["L0"][i]+sum(q[i,:,t]) <= data["L"][i]) #(25)
				@constraint(m,data["L0"][i]+sum(q[i,:,t])==data["d"][i-1,t]+I[i,t]) #22
			end
		end
		for k =1:data["k"]

			@constraint(m,sum(q[i,k,t] for i=2:data["n"]+1)<= data["Q"]*z[1,k,t]) #(30)
			@constraint(m,sum(x[1,j,k,t] for j=1:data["n"]+1 if j!=1)+sum(x[j,1,k,t] for j =1:data["n"]+1 if j!=1)==2*z[1,k,t]) #28
			
			for myi in [2:data["n"]+1;]
				@constraint(m,q[myi,k,t]<=Mit_til(data,myi-1,t)*z[myi,k,t]) #(26)
				@constraint(m,( sum(x[myi,j,k,t] for j =1:data["n"]+1 if j!=myi) + sum(x[j,myi,k,t] for j =1:data["n"]+1 if j!=myi) )== 2*z[myi,k,t]) #28
				
				#for s in subsets([2 : data["n"]+1],myi)
			end
			if(!branch)
				for s in powerset([ind for ind=2:data["n"]+1],2,data["n"])
					@constraint(m, sum(x[i,j,k,t] for i in s for j in s if i!=j) <= length(s)-1) #(29)
				end
			end
		
		end
	end
	return m
end

