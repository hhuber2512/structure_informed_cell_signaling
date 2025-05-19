using Distributions, DifferentialEquations
include("augment_likelihood.jl")
include("calculate_experimental_outputs.jl")
include("ode_problem.jl")
include("misc.jl")

function return_prior_bounds_empirical(bound_offset, binding_bounds_offset)
    values = groundtruth_parameter_values()
    prior_bounds = [[log10.(values[i])-bound_offset[1], log10.(values[i])+bound_offset[2]] for i in 1:50]
    koff_bounds = [-6-binding_bounds_offset[1],0+binding_bounds_offset[2]]
    kon_bounds = [-4-binding_bounds_offset[1],0+binding_bounds_offset[2]]
    kon_indices = return_kon_indices()
    koff_indices = return_koff_indices()
    for i in kon_indices
        prior_bounds[i] = kon_bounds
    end
    for i in koff_indices
        prior_bounds[i] = koff_bounds
    end
    return prior_bounds
end

"""

"""
function logprob_augmented(parameters, prior_dist, data, σ, index_training_data, ligand_dosages, odeproblems, odesolver_inputs, ode_solution_type, augment_likelihood)
    n_parameters = length(parameters)
    koff_indices = return_koff_indices()
    kon_indices = return_kon_indices()
    binding_reactions = return_order_of_binding_rxns()
    n_binding_reactions = length(binding_reactions)
    n_initial_conditions = length(odeproblems)

    #prior log probability
    logprior = 0.
    for i in 1:n_parameters
        logprior += logpdf(prior_dist[i], parameters[i])
    end

    #early exit if logprior = -Inf
    if logprior == -Inf
        logprob = -Inf
        return logprob
    end

    #log likelihood of augmented data 
    logaugmented = 0.
    for i in 1:n_binding_reactions
        kb = koff_indices[i]
        kf = kon_indices[i]
        logaugmented += regularize(augment_likelihood, parameters[kb], parameters[kf], binding_reactions[i])
    end

    #log likelihood of baseline data - must simulate ode for this
    #convert from log10 scale
    p_new = 10.0.^parameters
    #ode solver inputs
    save_at = odesolver_inputs["saveat"]
    solver = odesolver_inputs["solver"]
    abstol = odesolver_inputs["abstol"]
    reltol = odesolver_inputs["reltol"]
    #simulate for different initial conditions with same proposed parameter values
    predicted = Array{ode_sol_type}(undef,n_initial_conditions)
    for i in 1:n_initial_conditions
        op = remake(odeproblems[i], p=p_new)
        predicted[i] = DifferentialEquations.solve(op, solver, abstol=abstol, reltol=reltol, saveat=save_at)
    end
    #Early exit if simulation could not be computed successfully.
    if any([predicted[i].retcode !== ReturnCode.Success for i in 1:length(predicted)])
        logprob = -Inf
        return logprob
    end
    #calculate experimental quantities and reshape for likelihood
    experimental_quantities = [calculate_all_quantities(predicted[i]) for i in 1:n_initial_conditions]
    likelihood_quantities = [index_training_data(experimental_quantities[i], ligand_dosages[i]) for i in 1:n_initial_conditions]
    all_predictions_1d = Array{Float64}(undef,0)
    [append!(all_predictions_1d, likelihood_quantities[i]) for i in 1:n_initial_conditions]
    #data likelihood
    loglikelihood = logpdf(MvNormal(all_predictions_1d, σ), data) # MvNormal can take a vector input for standard deviation, and assumes covariance matrix is diagonal
    
    return logprior+logaugmented+loglikelihood
end

"""

"""
function logprob_baseline(parameters, prior_dist, data, σ, index_training_data, ligand_dosages, odeproblems, odesolver_inputs, ode_solution_type)
    n_parameters = length(parameters)
    koff_indices = return_koff_indices()
    kon_indices = return_kon_indices()
    binding_reactions = return_order_of_binding_rxns()
    n_binding_reactions = length(binding_reactions)
    n_initial_conditions = length(odeproblems)

    #prior log probability
    logprior = 0.
    for i in 1:n_parameters
        logprior += logpdf(prior_dist[i], parameters[i])
    end

    #early exit if logprior = -Inf
    if logprior == -Inf
        logprob = -Inf
        return logprob
    end

    #log likelihood of augmented data is set to zero in baseline case
    logaugmented = 0.

    #log likelihood of baseline data - must simulate ode for this
    #convert from log10 scale
    p_new = 10.0.^parameters
    #ode solver inputs
    save_at = odesolver_inputs["saveat"]
    solver = odesolver_inputs["solver"]
    abstol = odesolver_inputs["abstol"]
    reltol = odesolver_inputs["reltol"]
    #simulate for different initial conditions with same proposed parameter values
    predicted = Array{ode_sol_type}(undef,n_initial_conditions)
    for i in 1:n_initial_conditions
        op = remake(odeproblems[i], p=p_new)
        predicted[i] = DifferentialEquations.solve(op, solver, abstol=abstol, reltol=reltol, saveat=save_at)
    end
    #Early exit if simulation could not be computed successfully.
    if any([predicted[i].retcode !== ReturnCode.Success for i in 1:length(predicted)])
        logprob = -Inf
        return logprob
    end
    #calculate experimental quantities and reshape for likelihood
    experimental_quantities = [calculate_all_quantities(predicted[i]) for i in 1:n_initial_conditions]
    likelihood_quantities = [index_training_data(experimental_quantities[i], ligand_dosages[i]) for i in 1:n_initial_conditions]
    all_predictions_1d = Array{Float64}(undef,0)
    [append!(all_predictions_1d, likelihood_quantities[i]) for i in 1:n_initial_conditions]
    #data likelihood
    loglikelihood = logpdf(MvNormal(all_predictions_1d, σ), data) # MvNormal can take a vector input for standard deviation, and assumes covariance matrix is diagonal
    
    return logprior+logaugmented+loglikelihood
end