using JuMP
using CPLEX

function PL_LSP_Cvisite(n, l, f, u, d, h, L, M)
PB a partir de la formulation agrégée ??
    """Le PL du problème LSP en prenant en compte les coûts de visite
    Parametres :
    --- Notations LSP --- 
    n le nombre de revendeurs (en comptant le fournisseur)
    l horizon de planification
    f coût fixe par période de production
    u coût unitaire
    d tableau de dimension n*l, d[i][t] exprime la demande du revendeur i au temps t
    h tableau de dimension n, h[i] exprime le coût de stockage unitaire du revendeur i
    L tableau de dimension n, L[i] exprime la capacité de stockage maximale chez le revendeur i
    M tableau de dimension l, M[t] est une constante big M
    --- Notations VRP ---
    c tableau de dimension (n+1)*(n+1) est le coût de transport du sommet i vers le sommet j
    """

   # Création d'un modèle. Ce modèle fera l'interface avec le solveur GLPK
   m = Model(CPLEX.Optimizer)

   # Création des variables
   @variable(m, 0<=p[1:l])
   @variable(m, y[1:l], Bin)
   @variable(m, 0<=I[1:n][1:l])
   @variable(m, 0<=q[1:n][1:l])
   @variable(m, 0z[1:n][1:l], Bin)

   # Definition de Sc[i][t]
   SC = [c[0][i]+c[i][0] for i=1:n, t=1:l]
   
   # Fonction objectif
   @objective(m, Min, sum(u*p[t] + f*y[t] + sum(h[i]*I[i][t]+SC[i][t]*z[i][t] for i = 1:n) for t = 1:l ) )

   # Ajout des contraintes dans le modèle
   for t in 1:l
       @constraint(m, I[0][t-1] + p[t] = sum(q[i][t] for i = 1:n) + I[0][t] ) # Contrainte (1)
       @constraint(m, p[t] <= M[t]*y[t] # Contrainte (3)
       @constraint(m, I[0][t-1] <= L[0] # Contrainte (4)

       for i in 1:n
           @constraint(m, I[i][t-1] + q[i][t] = d[i][t] + I[i][t] ) # Contrainte (2)
           @constraint(m, I[i][t-1] + q[i][t] <= L[i] ) # Contrainte (5)

    # On n'en n'a pas besoin ici mais il peut être important de nommer les variables pour récupérer la solution duale
    #  Pour l'afficher, on fera println("\t c1 = ", -dual(c1)) # la fonction dual renvoie le coût réduit qui est l'opposé de la variable duale en maximisation

   end

# Affichages
   #Affichage du modèle
    println("Affichage du modèle avant résolution:")
    print(m)
    println()

    #Résolution du problème d'optimisation linéaire m par le solveur GLPK
    println("Résolution par le solveur linéaire choisi")
    optimize!(m)
    println()
   
    # Affiche tous les détails d'une solution à l'écran
    println("Affichage de tous les détails de la solution avec la commande solution_summary")
    println(solution_summary(m, verbose=true))
    println()

    # Mais on peut vouloir récupérer juste une information précise
    println("Récupération et affichage \"à la main\" d'informations précises")
    status = termination_status(m)

    if status == INFEASIBLE
        println("Le problème n'est pas réalisable")
    elseif status == UNBOUNDED
        println("Le problème est non borné")
    elseif status == OPTIMAL # ou JuMP.MathOptInterface.OPTIMAL
        println("Valeur optimale = ", objective_value(m))
        println("Solution optimale :")
        for i= 1:l
            println("\t p[",i,"] = ", value(p[i]))
        for i= 1:l
            println("\t y[",i,"] = ", value(y[i]))
        for i= 1:n
            for j= 1:l
                println("\t I[",i,"][",j,"] = ", value(I[i][j]))    
        for i= 1:n
            for j= 1:l
                println("\t q[",i,"][",j,"] = ", value(q[i][j]))
        for i= 1:n
            for j= 1:l
                println("\t SC[",i,"][",j,"] = ", SC[i][j])
        for i= 1:n
            for j= 1:l
                println("\t z[",i,"][",j,"] = ", value(z[i][j]))
        println("Temps de résolution :", solve_time(m))
    else
        println("Problème lors de la résolution")
    end
end
