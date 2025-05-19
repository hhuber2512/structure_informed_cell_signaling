using Plots, MCMCChains

function return_ess(chain, ndims, nwalkers, nsamples)
    #reshape from ndims, nwalkers, nsamples to nsamples, ndims, nwalkers
    #then, calculate ess
    ess = []
    for i in 1:ndims
        reshaped_chain = zeros(nsamples,1,nwalkers)
        for j in 1:nwalkers
            reshaped_chain[:,1,j] = chain[i,j,:]
        end
        mcmcchain = MCMCChains.Chains(reshaped_chain, ["param"])
        diagnostics = MCMCChains.ess(mcmcchain)
        ess = cat(ess, diagnostics[:,:ess], dims=1)
    end
    return ess
end

function trace_plots(rng, walkers_per_parameter, nwalkers, nsamples, parameter_name)
    n_walkers_to_plot = 5
    walkers_to_plot = Int.(ceil.(rand(rng, n_walkers_to_plot).*nwalkers)) #changes each time, so we get a look at more than the same 5 chains
    p = plot(1:nsamples, walkers_per_parameter[walkers_to_plot[1],:],label="walker $(walkers_to_plot[1])")
    for i in 2:n_walkers_to_plot
        plot!(1:nsamples, walkers_per_parameter[i,:],label="walker $(walkers_to_plot[i])")
    end
    xlabel!("iteration")
    ylabel!("$(parameter_name)")
    plot!(dpi=300, size=(690,690/1.5))
    return p
end
