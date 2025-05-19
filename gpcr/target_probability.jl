using Distributions, DifferentialEquations
include("augment_likelihood.jl")
include("calculate_experimental_outputs.jl")
include("ode_problem.jl")
include("misc.jl")

function return_prior_bounds_empirical(bound_offset, binding_bounds_offset)
    values = groundtruth_parameter_values()
    prior_bounds = [[log10.(values[i])-bound_offset[1], log10.(values[i])+bound_offset[2]] for i in 1:8]
    koff_bounds = [-6-binding_bounds_offset[1],0+binding_bounds_offset[2]]
    #must convert kon bounds from log10(1/(nM*sec)) to log10(1/(molecules*sec)) to match simulation units
    kon_bounds = [log10(convert_inverse_nM_molecules(10.0e-4))-binding_bounds_offset[1],log10(convert_inverse_nM_molecules(10.0e0))+binding_bounds_offset[2]]
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
function logprob_augmented(parameters, prior_dist, data, σ, odeproblems, odesolver_inputs_timecourse, odesolver_inputs_dose_response, ode_solution_type, augment_likelihood)
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
        #convert kf from log10(1/(molecules*sec)) to log10(1/(nM*sec)) for comparison 
        logaugmented += regularize(augment_likelihood, parameters[kb], log10(convert_inverse_molecules_nM(10.0^parameters[kf])), binding_reactions[i])
    end

    #log likelihood of baseline data - must simulate ode for this
    #convert from log10 scale
    p_new = 10.0.^parameters
    #ode solver inputs
    save_at_timecourse = odesolver_inputs_timecourse["saveat"]
    save_at_dose_response = odesolver_inputs_dose_response["saveat"]
    solver = odesolver_inputs_dose_response["solver"] #same for timecourse and dose response
    abstol = odesolver_inputs_dose_response["abstol"] #same for timecourse and dose response
    reltol = odesolver_inputs_dose_response["reltol"] #same for timecourse and dose response

    #first, simulate timecourse, using dose = 1000 nM (last ode problem) and save at = [0...600]
    predicted = DifferentialEquations.solve(remake(odeproblems[end], p=p_new), solver, abstol=abstol, reltol=reltol, saveat=save_at_timecourse);
    #early exit if simulation could not be computed successfully.
    if predicted.retcode !== ReturnCode.Success
        logprob = -Inf
        return logprob
    end
    #calculate experimental quantities
    active_g_timecourse = calculate_all_quantities(predicted)["fraction_active_G"]

    #next, simulate dose response
    predicted = Array{ode_sol_type}(undef, n_initial_conditions)
    for i in 1:n_initial_conditions
        op = remake(odeproblems[i], p=p_new)
        predicted[i] = DifferentialEquations.solve(op, solver, abstol=abstol, reltol=reltol, saveat=save_at_dose_response)
    end
    #early exit if simulation could not be computed successfully.
    if any([predicted[i].retcode !== ReturnCode.Success for i in 1:length(predicted)])
        logprob = -Inf
        return logprob
    end
    #calculate experimental quantities
    active_g = []
    for i in 1:n_initial_conditions
        fraction_activated = calculate_all_quantities(predicted[i])
        push!(active_g, fraction_activated["fraction_active_G"][1]) #need to index vector of 1d
    end
    active_g_dose_response = active_g[1:8]./active_g[end] #only have data for doses 0.1 - 100
    #reshape for likelihood
    all_predictions_1d = cat(active_g_timecourse, active_g_dose_response, dims=1)
    #data likelihood
    loglikelihood = logpdf(MvNormal(all_predictions_1d, σ), data) # MvNormal can take a vector input for standard deviation, and assumes covariance matrix is diagonal
    
    return logprior+logaugmented+loglikelihood
end

"""

"""
function logprob_baseline(parameters, prior_dist, data, σ, odeproblems, odesolver_inputs_timecourse, odesolver_inputs_dose_response, ode_solution_type)
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
    save_at_timecourse = odesolver_inputs_timecourse["saveat"]
    save_at_dose_response = odesolver_inputs_dose_response["saveat"]
    solver = odesolver_inputs_dose_response["solver"] #same for timecourse and dose response
    abstol = odesolver_inputs_dose_response["abstol"] #same for timecourse and dose response
    reltol = odesolver_inputs_dose_response["reltol"] #same for timecourse and dose response

    #first, simulate timecourse, using dose = 1000 nM (last ode problem) and save at = [0...600]
    predicted = DifferentialEquations.solve(remake(odeproblems[end], p=p_new), solver, abstol=abstol, reltol=reltol, saveat=save_at_timecourse);
    #early exit if simulation could not be computed successfully.
    if predicted.retcode !== ReturnCode.Success
        logprob = -Inf
        return logprob
    end
    #calculate experimental quantities
    active_g_timecourse = calculate_all_quantities(predicted)["fraction_active_G"]

    #next, simulate dose response
    predicted = Array{ode_sol_type}(undef, n_initial_conditions)
    for i in 1:n_initial_conditions
        op = remake(odeproblems[i], p=p_new)
        predicted[i] = DifferentialEquations.solve(op, solver, abstol=abstol, reltol=reltol, saveat=save_at_dose_response)
    end
    #early exit if simulation could not be computed successfully.
    if any([predicted[i].retcode !== ReturnCode.Success for i in 1:length(predicted)])
        logprob = -Inf
        return logprob
    end
    #calculate experimental quantities
    active_g = []
    for i in 1:n_initial_conditions
        fraction_activated = calculate_all_quantities(predicted[i])
        push!(active_g, fraction_activated["fraction_active_G"][1]) #need to index vector of 1d
    end
    active_g_dose_response = active_g[1:8]./active_g[end] #only have data for doses 0.1 - 100
    #reshape for likelihood
    all_predictions_1d = cat(active_g_timecourse, active_g_dose_response, dims=1)
    #data likelihood
    loglikelihood = logpdf(MvNormal(all_predictions_1d, σ), data) # MvNormal can take a vector input for standard deviation, and assumes covariance matrix is diagonal
    
    return logprior+logaugmented+loglikelihood
end

