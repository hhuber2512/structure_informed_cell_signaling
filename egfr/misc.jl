function return_kon_indices()
    return [10,14,16,22,26,28,32,46,52].-2
end

function return_koff_indices()
    return [9,13,15,21,25,27,31,45,51].-2
end

function return_order_of_binding_rxns()
    return ["k10", "k12", "k13", "k17", "k19", "k1", "k21", "k5", "k9"]
end

function return_nonbinding_indices()
    n_parameters = 50
    kon_indices = return_kon_indices()
    koff_indices = return_koff_indices()
    binding_indices = cat(return_koff_indices(), return_kon_indices(), dims=1)
    n_binding_parameters = length(binding_indices)
    return collect(1:n_parameters)[[all(fill(collect(1:n_parameters)[i],n_binding_parameters) .!= binding_indices) for i in 1:n_parameters]]
end

function species_indices_map()
    species_ordered = ["EGF", "R", "Ra", "R2", "RP", "RPLCg", "RPLCgP", "RG", "RGS", "RSh", "RShP", "RShG", "RShGS", "GS", "ShP","ShG","PLCg","PLCgP", "PLCgl","Grb","Shc","SOS","ShGS"]
    indices = collect(1:23)
    return Dict(species_ordered[i]=>indices[i] for i in 1:23)
end

function indices_species_map()
    species_ordered = ["EGF", "R", "Ra", "R2", "RP", "RPLCg", "RPLCgP", "RG", "RGS", "RSh", "RShP", "RShG", "RShGS", "GS", "ShP","ShG","PLCg","PLCgP", "PLCgl","Grb","Shc","SOS","ShGS"]
    indices = collect(1:23)
    return Dict(indices[i]=>species_ordered[i] for i in 1:23)
end

function parameters_indices_map()
    parameters_ordered = ["K16"
    "K4"
    "K8"
    "V16"
    "V4"
    "V8"
    "k10b"
    "k10f"
    "k11b"
    "k11f"
    "k12b"
    "k12f"
    "k13b"
    "k13f"
    "k14b"
    "k14f"
    "k15b"
    "k15f"
    "k17b"
    "k17f"
    "k18b"
    "k18f"
    "k19b"
    "k19f"
    "k1b"
    "k1f"
    "k20b"
    "k20f"
    "k21b"
    "k21f"
    "k22b"
    "k22f"
    "k23b"
    "k23f"
    "k24b"
    "k24f"
    "k25b"
    "k25f"
    "k2b"
    "k2f"
    "k3b"
    "k3f"
    "k5b"
    "k5f"
    "k6b"
    "k6f"
    "k7b"
    "k7f"
    "k9b"
    "k9f"]
    indices = collect(1:50)
    return Dict(parameters_ordered[i]=>indices[i] for i in 1:50)
end

function return_inferred_parameters()
    parameters_ordered = ["K16"
    "K4"
    "K8"
    "V16"
    "V4"
    "V8"
    "k10b"
    "k10f"
    "k11b"
    "k11f"
    "k12b"
    "k12f"
    "k13b"
    "k13f"
    "k14b"
    "k14f"
    "k15b"
    "k15f"
    "k17b"
    "k17f"
    "k18b"
    "k18f"
    "k19b"
    "k19f"
    "k1b"
    "k1f"
    "k20b"
    "k20f"
    "k21b"
    "k21f"
    "k22b"
    "k22f"
    "k23b"
    "k23f"
    "k24b"
    "k24f"
    "k25b"
    "k25f"
    "k2b"
    "k2f"
    "k3b"
    "k3f"
    "k5b"
    "k5f"
    "k6b"
    "k6f"
    "k7b"
    "k7f"
    "k9b"
    "k9f"]
    return parameters_ordered
end