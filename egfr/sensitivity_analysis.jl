using DifferentialEquations
include("target_probability.jl")


function local_sensitivity_analysis(p, prob, solver_inputs,species_index)
    p_new = 10.0.^p
    #p_new = p
    prob = remake(prob, p = p_new)
    solve(prob, solver_inputs["solver"], reltol=solver_inputs["reltol"], 
    abstol=solver_inputs["abstol"], saveat = solver_inputs["saveat"])[species_index,:]
end

function local_sensitivity_analysis_kd(p, prob, solver_inputs, species_index, kf, kb_indices, kf_indices, nonbinding_indices)
    n_kd = length(kf)
    #convert back to kb and kf form
    kb = p[1:n_kd].*kf 
    p_new = zeros(50)
    p_new[kb_indices] = kb
    p_new[kf_indices] = kf
    p_new[nonbinding_indices] = p[n_kd+1:end]
    p_new = 10.0.^p_new
    prob = remake(prob, p = p_new)
    solve(prob, solver_inputs["solver"], reltol=solver_inputs["reltol"], 
    abstol=solver_inputs["abstol"], saveat = solver_inputs["saveat"])[species_index,:]
end

function global_sensitivity_analysis(p, ligand_doses, prob, solver_inputs, n_species, non_zero_ic_mask)
    p = 10.0.^p
    u0_new = p[1:5] #initial conditions
    p_new = p[6:end] #parameters
    all_sim = [0.0]
    for j in 1:length(ligand_doses)
        dummy1, u0, dummy2, dummy3 = return_ode_problem_default_inputs(ligand_doses[j])
        u0[non_zero_ic_mask] = u0_new 
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

function global_sensitivity_analysis_mi(p, ligand_dose, prob, solver_inputs, species_index, time_index)
    p_new = 10.0.^p
    dummy1, u0, dummy2, dummy3 = return_ode_problem_default_inputs(ligand_dose)
    prob = remake(prob, p = p_new, u0 = u0)
    output = solve(prob, solver_inputs["solver"], reltol=solver_inputs["reltol"], 
    abstol=solver_inputs["abstol"], saveat = solver_inputs["saveat"])
    return output[species_index, time_index]
end

function sensitivity_analysis_bounds(bound_offset,binding_bound_offset)
    prior_bounds = return_prior_bounds_empirical(bound_offset, binding_bound_offset) #for target probability
    dummy = [-4,0]
    return push!(prior_bounds, dummy)
end

function save_global_sensitivity_indices_1d(path)
    result = deserialize("$(path).jls") #exclude initial conditions
    #save sum of S1 w.r.t. different parameter types
    global_ic = sum(result.S1[:,1:5], dims=2)
    global_kb = sum(result.S1[:,return_koff_indices().+5], dims=2) #shift to account for sensitivity of initial values, which we exclude from our analysis
    global_kf = sum(result.S1[:,return_kon_indices().+5], dims=2)
    global_nonbind = sum(result.S1[:,return_nonbinding_indices().+5], dims=2)
    global_dummy = sum(result.S1[:,51+5], dims=2) #dummy variable sensitivity
    serialize("$(path)_S1_ic.jls", global_ic) 
    serialize("$(path)_S1_kb.jls", global_kb)
    serialize("$(path)_S1_kf.jls", global_kf)
    serialize("$(path)_S1_nonbind.jls", global_nonbind)
    serialize("$(path)_S1_dummy.jls", global_dummy)
    #save max of ST w.r.t. different parameter types
    global_ic = maximum(result.ST[:,1:5], dims=2)
    global_kb = maximum(result.ST[:,return_koff_indices().+5], dims=2)
    global_kf = maximum(result.ST[:,return_kon_indices().+5], dims=2)
    global_nonbind = maximum(result.ST[:,return_nonbinding_indices().+5], dims=2)
    global_dummy = maximum(result.ST[:,51+5], dims=2) #dummy variable sensitivity 
    serialize("$(path)_ST_ic.jls", global_ic) 
    serialize("$(path)_ST_kb.jls", global_kb)
    serialize("$(path)_ST_kf.jls", global_kf)
    serialize("$(path)_ST_nonbind.jls",global_nonbind) 
    serialize("$(path)_ST_dummy.jls", global_dummy)
end
