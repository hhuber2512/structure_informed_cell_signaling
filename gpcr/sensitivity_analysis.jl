using DifferentialEquations


function local_sensitivity_analysis(p, prob, solver_inputs,species_index)
    p_new = 10.0.^p
    prob = remake(prob, p = p_new)
    solve(prob, solver_inputs["solver"], reltol=solver_inputs["reltol"], 
    abstol=solver_inputs["abstol"], saveat = solver_inputs["saveat"])[species_index,:]
end

function global_sensitivity_analysis(p, prob, solver_inputs, n_species)
    p_new = 10.0.^p
    prob = remake(prob, p = p_new)
    output = solve(prob, solver_inputs["solver"], reltol=solver_inputs["reltol"], 
    abstol=solver_inputs["abstol"], saveat = solver_inputs["saveat"])
    flattened_simulation = [0.0] #initialize as Float64
    for i in 1:n_species
        flattened_simulation = cat(flattened_simulation, output[i,:], dims=1)
    end
    return flattened_simulation[2:end]
end