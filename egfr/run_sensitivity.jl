using Serialization, Statistics, Distributions, StableRNGs, CSV, ForwardDiff
include("misc.jl")
include("ode_problem.jl")
include("target_probability.jl")
include("marginal_kl_divergence.jl")
include("define_likelihood_data.jl")
include("sensitivity_analysis.jl");

#for analysis needing pairs (augmented and baseline comparison)
augmented_run = "220"
baseline_run = "201"
prior_nonbinding_offset = [2,2]
prior_binding_offset = [0,0];

#first, extract parameter set that maximizes baseline training data likelihood
likelihood = deserialize("outputs/300_$(baseline_run)_posterior_samples_thinned_likelihood.jls");
max_likelihood = maximum(likelihood)
mask = likelihood .== max_likelihood
max_likelihood_parameters = deserialize("outputs/300_$(baseline_run)_posterior_samples_thinned.jls")[:, mask]

#calculate local sensitivity of each species about this parameter value 
egf_doses = return_ligand_dose_order_for_likelihood_w_pshc_test()
n_doses = length(egf_doses)
n_species = 23
solver_inputs = return_ode_problem_solver_default_inputs()
#redefined save at to be finer grain
#solver_inputs["saveat"] = return_finegrain_saveat()
local_sensitivity_kb = []
local_sensitivity_kf = []
local_sensitivity_nonbinding = []
for i in 1:n_doses
    odesys, u0, tspan, p = return_ode_problem_default_inputs(egf_doses[i])
    odeprob = DifferentialEquations.ODEProblem(odesys, u0, tspan, p) #note, parameters will be redefined in next function
    for j in 1:n_species
        local_sensitivity = ForwardDiff.jacobian(p -> local_sensitivity_analysis(p, odeprob, solver_inputs, j), max_likelihood_parameters)
        push!(local_sensitivity_kb, [local_sensitivity[:,return_koff_indices()]])
        push!(local_sensitivity_kf, [local_sensitivity[:, return_kon_indices()]])
        push!(local_sensitivity_nonbinding, [local_sensitivity[:, return_nonbinding_indices()]])
    end
end
serialize("outputs/600_$(baseline_run)_local_sensitivity_kb.jls", local_sensitivity_kb)
serialize("outputs/600_$(baseline_run)_local_sensitivity_kf.jls", local_sensitivity_kf)
serialize("outputs/600_$(baseline_run)_local_sensitivity_nonbind.jls", local_sensitivity_nonbinding)
