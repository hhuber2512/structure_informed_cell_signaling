using Serialization, AffineInvariantMCMC

function combine_and_thin_posterior_samples(my_run, thin_each_walker_by, n_subchains)
    #exclude first 1000 chains as burn-in
    #thin first and second half for memory purposes
    chain = deserialize("outputs/$(my_run)_posterior_samples_subchain_1.jls")
    for i in 2:Int(n_subchains/2)
        subchain = deserialize("outputs/$(my_run)_posterior_samples_subchain_$(i).jls")
        chain = cat(chain, subchain, dims=3)
    end
    chain = chain[:,:,1:thin_each_walker_by:end]; #ndims x nwalkers x nsamplers_per_walker
    
    chain2 = deserialize("outputs/$(my_run)_posterior_samples_subchain_$((Int(n_subchains/2)+1)).jls")
    for i in (Int(n_subchains/2)+2):n_subchains
        subchain = deserialize("outputs/$(my_run)_posterior_samples_subchain_$(i).jls")
        chain2 = cat(chain2, subchain, dims=3)
    end

    chain2 = chain2[:,:,1:thin_each_walker_by:end]; #ndims x nwalkers x nsamplers_per_walker

    chain = cat(chain, chain2, dims=3)

    ll = deserialize("outputs/$(my_run)_posterior_samples_likelihood_subchain_1.jls")
    for i in 2:n_subchains
        sub_ll = deserialize("outputs/$(my_run)_posterior_samples_likelihood_subchain_$(i).jls")
        ll = cat(ll,sub_ll, dims=2)
    end
    ll = ll[:,1:thin_each_walker_by:end]; #nwalkers x nsamples
    return AffineInvariantMCMC.flattenmcmcarray(chain, ll)
end

function combine_posterior_samples_keep_walkers(my_run, thin_each_walker_by, n_subchains, p_index)
    #exclude first 1000 chains as burn-in
    #take portion of parmeters per p_index
    chain = deserialize("outputs/$(my_run)_posterior_samples_subchain_1.jls")[p_index,:,:]
    for i in 2:n_subchains
        subchain = deserialize("outputs/$(my_run)_posterior_samples_subchain_$(i).jls")[p_index,:,:]
        chain = cat(chain, subchain, dims=3)
    end
    chain1 = chain[:,:,1:thin_each_walker_by:end]; #ndims x nwalkers x nsamplers_per_walker

    ll = deserialize("outputs/$(my_run)_posterior_samples_likelihood_subchain_1.jls")
    for i in 2:n_subchains
        sub_ll = deserialize("outputs/$(my_run)_posterior_samples_likelihood_subchain_$(i).jls")
        ll = cat(ll,sub_ll, dims=2)
    end
    ll = ll[:,1:thin_each_walker_by:end]; #nwalkers x nsamples
    return chain, ll
end