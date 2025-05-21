using DifferentialEquations


function local_sensitivity_analysis(p, prob, solver_inputs,species_index)
    p_new = 10.0.^p
    prob = remake(prob, p = p_new)
    solve(prob, solver_inputs["solver"], reltol=solver_inputs["reltol"], 
    abstol=solver_inputs["abstol"], saveat = solver_inputs["saveat"])[species_index,:]
end

#function global_sensitivity_analysis(p, prob, solver_inputs, n_species)
    #p_new = 10.0.^p
    #prob = remake(prob, p = p_new)
    #output = solve(prob, solver_inputs["solver"], reltol=solver_inputs["reltol"], 
    #abstol=solver_inputs["abstol"], saveat = solver_inputs["saveat"])
    #flattened_simulation = [0.0] #initialize as Float64
    #for i in 1:n_species
        #flattened_simulation = cat(flattened_simulation, output[i,:], dims=1)
    #end
    #return flattened_simulation[2:end]
#end

function global_sensitivity_analysis(p, ligand_doses, prob, solver_inputs, n_species)
    p_new = 10.0.^p
    all_sim = [0.0]
    for j in 1:length(ligand_doses)
        dummy1, u0, dummy2, dummy3 = return_ode_problem_default_inputs(ligand_doses[j])
        prob = remake(prob, p = p_new, u0 = u0)
        output = solve(prob, solver_inputs["solver"], reltol=solver_inputs["reltol"], 
        abstol=solver_inputs["abstol"], saveat = solver_inputs["saveat"])
        flattened_simulation = [0.0] #initialize as Float64
        for i in 1:n_species
            flattened_simulation = cat(flattened_simulation, output[i,:], dims=1)
        end
        all_sim = cat(all_sim, flattened_simulation[2:end], dims=1) #leave out initialization of array flattened_sim...
    end
    return all_sim[2:end]
end

function save_global_sensitivity_indices_1d(path)
    result = deserialize("$(path).jls")
    #save sum of S1 w.r.t. different parameter types
    global_kb = sum(result.S1[:,return_koff_indices()], dims=2)
    global_kf = sum(result.S1[:,return_kon_indices()], dims=2)
    global_nonbind = sum(result.S1[:,return_nonbinding_indices()], dims=2)
    global_dummy = sum(result.S1[:,9], dims=2) #dummy variable sensitivity 
    serialize("$(path)_S1_kb.jls", global_kb)
    serialize("$(path)_S1_kf.jls", global_kf)
    serialize("$(path)_S1_nonbind.jls", global_nonbind)
    serialize("$(path)_S1_dummy.jls", global_dummy)
    #save max of ST w.r.t. different parameter types
    global_kb = maximum(result.ST[:,return_koff_indices()], dims=2)
    global_kf = maximum(result.ST[:,return_kon_indices()], dims=2)
    global_nonbind = maximum(result.ST[:,return_nonbinding_indices()], dims=2)
    global_dummy = maximum(result.ST[:,9], dims=2) #dummy variable sensitivity 
    serialize("$(path)_ST_kb.jls", global_kb)
    serialize("$(path)_ST_kf.jls", global_kf)
    serialize("$(path)_ST_nonbind.jls",global_nonbind) 
    serialize("$(path)_ST_dummy.jls", global_dummy)
end