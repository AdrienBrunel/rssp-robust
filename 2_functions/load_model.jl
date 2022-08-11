# ==============================================================================
#  NOMINAL RESERVE SITE SELECTION ; MINSET ; WITH BLM ; WITH LINEARIZATION
# ==============================================================================
	function reserve_site_selection_model(instance,gridgraph)

		# Lecture données
		Arcs                 = gridgraph.Arcs
		NoeudsPeripheriques  = gridgraph.NoeudsPeripheriques
		ConservationFeatures = instance.ConservationFeatures
		PlanningUnits 		 = instance.PlanningUnits
		LockedOut 		     = instance.LockedOut
		Amount 				 = instance.Amount
		Cost                 = instance.Cost
		Targets              = instance.Targets
		BoundaryLength       = instance.BoundaryLength
		BoundaryCorrection   = instance.BoundaryCorrection
		Beta   				 = instance.Beta

		# Déclaration du modèle
		m = Model(Cbc.Optimizer)
		set_optimizer_attribute(m, "logLevel", 1)
		set_optimizer_attribute(m, "seconds", 600)
		#m = Model(Gurobi.Optimizer)
		#set_optimizer_attribute(m, "TimeLimit", 600)

		# Variables de décision
		@variable(m, x[PlanningUnits], Bin); # variable de sélection du noeuds j dans le graphe de la réserve
		@variable(m, z[Arcs], Bin);          # variable de linéarisation du périmètre quadratique

		# Minimiser le coût de la réserve
		@objective(m, Min, sum(Cost[j]*x[j] for j in PlanningUnits) + Beta*sum(BoundaryLength[d]*(x[d[1]]-z[d]) for d in Arcs) + Beta*sum(BoundaryCorrection[j]*x[j] for j in NoeudsPeripheriques))

		# La réserve doit réaliser les cibles écologiques
		@constraint(m, cibles[i in ConservationFeatures], sum(Amount[i,j]*x[j] for j in PlanningUnits) >= Targets[i])

		# linearization constraints
		@constraint(m, lin_z1[d in Arcs],  z[d] - x[d[1]] <= 0)
		@constraint(m, lin_z2[d in Arcs],  z[d] - x[d[2]] <= 0)

		# locked out planning units constraints
		@constraint(m, locked_out[j in LockedOut],  x[j] == 0)

		return m
	end


# ==============================================================================
#  ROBUST RESERVE SITE SELECTION ; MINSET ; WITH BLM ; WITH LINEARIZATION 
# ==============================================================================
	function minset_robust(instance,gridgraph,Gamma)

		# Loading data
		Arcs                 = gridgraph.Arcs
		NoeudsPeripheriques  = gridgraph.NoeudsPeripheriques
		ConservationFeatures = instance.ConservationFeatures
		PlanningUnits 		 = instance.PlanningUnits
		LockedOut 		     = instance.LockedOut
		Amount 				 = instance.Amount
		Sigma 				 = instance.Sigma
		Cost                 = instance.Cost
		Targets              = instance.Targets
		BoundaryLength       = instance.BoundaryLength
		BoundaryCorrection   = instance.BoundaryCorrection
		Beta   				 = instance.Beta

		# Model declaration
		m = Model(Cbc.Optimizer)
		set_optimizer_attribute(m, "logLevel", 1)
		set_optimizer_attribute(m, "seconds", 600)
		#set_optimizer_attribute(m, "maxnodes", 100000)
		#set_optimizer_attribute(m, "ratiogap", 0.001)

		# Decision variables
		@variable(m, x[PlanningUnits], Bin); # 1 if planning unit j in reserve 
		@variable(m, z[Arcs], Bin);          # linearization variable for x_j*x_i
		@variable(m, u[ConservationFeatures,PlanningUnits] >= 0) # robust optimization dual variable
		@variable(m, v[ConservationFeatures] >= 0) # robust optimization dual variable

		# Minset objective
		@objective(m, Min, sum(Cost[j]*x[j] for j in PlanningUnits) + Beta*sum(BoundaryLength[d]*(x[d[1]]-z[d]) for d in Arcs) + Beta*sum(BoundaryCorrection[j]*x[j] for j in NoeudsPeripheriques))

		# The reserve must meet the robust conservation features targets 
		@constraint(m, targets[i in ConservationFeatures],  sum(Amount[i,j]*x[j] for j in PlanningUnits) - (Gamma[i]*v[i]+sum(u[i,j] for j in PlanningUnits)) >= Targets[i])

		# locked out/in constraints
		@constraint(m, locked_out[j in LockedOut],  x[j] == 0)

		# linearization constraints
		@constraint(m, lin_z1[d in Arcs],  z[d] - x[d[1]] <= 0)
		@constraint(m, lin_z2[d in Arcs],  z[d] - x[d[2]] <= 0)

		# robust constraint
		@constraint(m, rob_uv[i in ConservationFeatures,j in PlanningUnits],  v[i] + u[i,j] >= x[j]*Sigma[i,j])

		return m
	end
