using Serialization, AffineInvariantMCMC

function combine_and_thin_posterior_samples(my_run, thin_each_walker_by, n_subchains)
    #exclude first 1000 chains as burn-in
    chain = deserialize("outputs/$(my_run)_posterior_samples_subchain_1.jls")
    for i in 2:n_subchains
        subchain = deserialize("outputs/$(my_run)_posterior_samples_subchain_$(i).jls")
        chain = cat(chain, subchain, dims=3)
    end
    chain = chain[:,:,1:thin_each_walker_by:end]; #ndims x nwalkers x nsamplers_per_walker

    ll = deserialize("outputs/$(my_run)_posterior_samples_likelihood_subchain_1.jls")
    for i in 2:n_subchains
        sub_ll = deserialize("outputs/$(my_run)_posterior_samples_likelihood_subchain_$(i).jls")
        ll = cat(ll,sub_ll, dims=2)
    end
    ll = ll[:,1:thin_each_walker_by:end]; #nwalkers x nsamples
    return AffineInvariantMCMC.flattenmcmcarray(chain, ll)
end

function combine_and_thin_posterior_samples_keep_walkers(my_run, thin_each_walker_by, n_subchains)
    #exclude first 1000 chains as burn-in
    chain = deserialize("outputs/$(my_run)_posterior_samples_subchain_1.jls")
    for i in 2:n_subchains
        subchain = deserialize("outputs/$(my_run)_posterior_samples_subchain_$(i).jls")
        chain = cat(chain, subchain, dims=3)
    end
    chain = chain[:,:,1:thin_each_walker_by:end]; #ndims x nwalkers x nsamplers_per_walker

    ll = deserialize("outputs/$(my_run)_posterior_samples_likelihood_subchain_1.jls")
    for i in 2:n_subchains
        sub_ll = deserialize("outputs/$(my_run)_posterior_samples_likelihood_subchain_$(i).jls")
        ll = cat(ll,sub_ll, dims=2)
    end
    ll = ll[:,1:thin_each_walker_by:end]; #nwalkers x nsamples
    return chain, ll
end