using IterTools
using JuMP
using CPLEX
include("PDI1.jl")
include("PDI2.jl")
include("toolsPDI.jl")
include("tools.jl")


#---------------------------------------------------------------
#PRP 2
#---------------------------------------------------------------
function BC2(data)
	"""
	data : donnees des instances sous forme de dictionnaire
	renvoie la solution optimale , utilisation du branch and cut et du PLNE PDI2
	"""

	#m=PLNE_PDI2(data,true) #model du PLNE PDI2
	branch=true
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
			for i =1:data["n"]+1
				@constraint(m,x[i,i,k,t]==0)
			end
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
	
	function lazySep_Violated(cb_data)
		for t = 1:data["l"]
			for k =1:data["k"]
				xsep=zeros(Int,data["n"]+1,data["n"]+1)
				for i=1:data["n"]+1
					for j =1:data["n"]+1
						if( callback_value(cb_data,x[i,j,k,t]) > 0.9)
							xsep[i,j]=1
						end
					end
				end
				println(xsep)
				departDepot,tournee=detect_sous_tour2(data["n"],xsep)
				println(departDepot,tournee)
				
				if(!departDepot)
					if(length(tournee)>=2 && !(1 in tournee))
						for ee in tournee
							#con=@build_constraint(sum(variable_by_name(m,"x[$i,$j,$k,$t]") for i in tournee for j in tournee if i!=j ) <= sum(variable_by_name(m,"z[$i,$k,$t]") for i in tournee) -variable_by_name(m,"z[$e,$k,$t]"))
							con=@build_constraint(sum(x[i,j,k,t] for i in tournee for j in tournee if i!=j ) <= sum(z[i-1,k,t] for i in tournee) -z[ee-1,k,t])
							MOI.submit(m,MOI.LazyConstraint(cb_data),con)
						end
					end
				end
			end
			
		end
	end
	
	
	MOI.set(m,MOI.LazyConstraintCallback(),lazySep_Violated)
	
	optimize!(m)
	
	status = termination_status(m)
	if status == JuMP.MathOptInterface.INFEASIBLE
		println("Le problème n'est pas réalisable")
	elseif status == JuMP.MathOptInterface.DUAL_INFEASIBLE
		println("Le problème est non borné")
	elseif status == JuMP.MathOptInterface.OPTIMAL
		xopt,yopt,qopt,popt=get_param_opt2(data,m)
		println("Temps de résolution :", solve_time(m))
		println("optimum = ", objective_value(m)) 
		return true,objective_value(m), solve_time(m) ,xopt,yopt,qopt,popt,relative_gap(m)
	end
	false,0,0,0,0,0,0,0

end



#---------------------------------------------------------------
#PRP 1
#---------------------------------------------------------------

function BC1(data)
	"""
	data : donnees des instances sous forme de dictionnaire
	renvoie la solution optimale , utilisation du branch and cut et du PLNE PDI1
	"""
	#m=PLNE_PDI1(data,true)
	branch=true
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
	@variable(m,0 <= q[1:data["n"]+1,1:data["l"]],Int) #quantite delivree au revendeur i a la periode tincl
	@variable(m, x[1:data["n"]+1,1:data["n"]+1,1:data["l"]],Bin) #x_ijt=1 si un vehicule voyage directement du noeud i au noeud j a l instant t
	
	@objective(m, Min, sum(data["u"]*p + data["f"]*y + sum( data["h"] .*I,dims=1 )[1,:] + sum(sum(lci .*x, dims=1), dims=2)[1,1,:] )) 
	
	#contraintes
	for t in [1:data["l"];]
		for i =1:data["n"]+1
			@constraint(m,x[i,i,t]==0)
		end
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
	
	function lazySep_Violated(cb_data)
		for t = 1:data["l"]
			xsep=zeros((data["n"]+1,data["n"]+1))
			for i=1:data["n"]+1
				for j =1:data["n"]+1
					if( round(callback_value(cb_data,x[i,j,t]))==1)
						xsep[i,j]=1
					end
				end
			end
			trouveDepartDepot,sousTours=detect_sous_tour1(data["n"],xsep)
			if(!trouveDepartDepot) # trouver au moins un sous tours
				for tournee in sousTours
					if(length(tournee)>=2)
						#con = @build_constraint(sum(variable_by_name(m,"x[$i,$j,$t]") for i in tournee for j in tournee if i!=j )*data["Q"] <= sum( data["Q"]*variable_by_name(m,"z[$i,$t]") - variable_by_name(m,"q[$i,$t]") for i in tournee))
						con = @build_constraint(sum(x[i,j,t] for i in tournee for j in tournee if i!=j )*data["Q"] <= sum( data["Q"]*z[i-1,t] - q[i,t] for i in tournee))
						MOI.submit(m,MOI.LazyConstraint(cb_data),con)
						sbar=setdiff(Set(1:data["n"]+1),tournee)
						con = @build_constraint(sum(x[i,j,t] for i in sbar for j in tournee) >= sum(q[i,t] for i in tournee)/data["Q"]) 
						MOI.submit(m,MOI.LazyConstraint(cb_data),con)
			
					end
				end
			end
		end
	end
	set_time_limit_sec(m,600)
	MOI.set(m,MOI.LazyConstraintCallback(),lazySep_Violated)
	optimize!(m)
	
	status = termination_status(m)
	if status == JuMP.MathOptInterface.INFEASIBLE
		println("Le problème n'est pas réalisable")
	elseif status == JuMP.MathOptInterface.DUAL_INFEASIBLE
		println("Le problème est non borné")
	elseif status == JuMP.MathOptInterface.OPTIMAL
		xopt,yopt,qopt,popt=get_param_opt1(data,m)
		println("Temps de résolution :", solve_time(m))
		println("optimum = ", objective_value(m)) 
		return true,objective_value(m), solve_time(m) ,xopt,yopt,qopt,popt,relative_gap(m)
		#println(solution_summary(m, verbose=true))
		#return true,
	end
	return false,0,0,0,0,0,0,0
	

end



