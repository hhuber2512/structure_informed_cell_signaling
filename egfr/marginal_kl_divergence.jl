using KernelDensity, Distributions

function kl_divergence_bits_samples(p,q,parameter_names)
    #take in arrays, calculate kl divergence
    #how much information do we lose using q to approximate p
    #approximate kl divergence as expectations w.r.t. p
    n_parameters, n_iterations = size(p)

    kl_div = Array{Float64}(undef, n_parameters)
    for i in 1:n_parameters
        U_q = kde(q[i,:])
        ik_q = InterpKDE(U_q)
        U_p = kde(p[i,:])
        ik_p = InterpKDE(U_p)
        num = pdf(ik_p, p[i,:])
        denom = pdf(ik_q, p[i,:]) #approximate as expectations w.r.t. p, take pdf of p samples
        #check for zeros, negative values, input very small number 1.0e-20 if true
        check = denom .<= 0
        denom[check] .= 1.0e-20
        check = num .<= 0
        num[check] .= 1.0e-20
        kl_div[i] = (1/n_iterations)*sum(log2.(num) .- log2.(denom))
    end
    return Dict(parameter_names[i]=> kl_div[i] for i in 1:n_parameters), parameter_names
end
