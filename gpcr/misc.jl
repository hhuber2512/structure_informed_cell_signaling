function return_kon_indices()
    return [1]
end

function return_koff_indices()
    return [2]
end

function return_order_of_binding_rxns()
    return ["k1"]
end

function return_nonbinding_indices()
    return [3,4,5,6,7,8]
end

function species_indices_map()
    species_ordered = ["R", "L", "RL", "G", "Ga"]
    indices = collect(1:5)
    return Dict(species_ordered[i]=>indices[i] for i in 1:5)
end

function indices_species_map()
    species_ordered = ["R", "L", "RL", "G", "Ga"]
    indices = collect(1:5)
    return Dict(indices[i]=>species_ordered[i] for i in 1:5)
end

function parameters_indices_map()
    parameters_ordered = ["k1f", "k1b", "k2", "k3", "k4", "k5", "k6", "k7"]
    indices = collect(1:8)
    return Dict(parameters_ordered[i]=>indices[i] for i in 1:8)
end

function return_inferred_parameters()
    parameters_ordered = ["k1f", "k1b", "k2", "k3", "k4", "k5", "k6", "k7"]
    return parameters_ordered
end

function convert_molecules_nM(x_molecules)
    #N.B. converting from molecules per cell to nM
    cell_volume = 8.6e-14 #Liters, average cell volume from bionumbers - https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=4&id=111978
    avo_constant = 6.022e23 #avogadros constant, molecules/molecules
    M_to_nM = 1e9 #nano conversion,  10^9
    nM = (x_molecules*M_to_nM)/(cell_volume*avo_constant)
    return nM
end

function convert_nM_molecules(x_nanomolar)
    #N.B. converting from nM to molecules per cell
    cell_volume = 8.6e-14 #Liters, average yeast cell volume from bionumbers - https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=4&id=111978
    avo_constant = 6.022e23 #avogadros constant, molecules/molecules
    M_to_nM = 1e9 #nano conversion, 10^9
    molecules_per_cell = (x_nanomolar*avo_constant*cell_volume)/M_to_nM
    return molecules_per_cell
end

function convert_inverse_molecules_nM(x_molecules)
    #N.B. converting from 1/(molecules per cell * sec) to 1/(nM*sec)
    cell_volume = 8.6e-14 #Liters, average cell volume from bionumbers - https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=4&id=111978
    avo_constant = 6.022e23 #avogadros constant, molecules/molecules
    M_to_nM = 1e9 #nano conversion, 10^9
    nM = (x_molecules*avo_constant*cell_volume)/(M_to_nM)
    return nM
end

function convert_inverse_nM_molecules(x_nanomolar)
    #N.B. converting from 1/(nM*sec) to 1/(molecules per cell * sec)
    cell_volume = 8.6e-14 #Liters, average cell volume from bionumbers - https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=4&id=111978
    avo_constant = 6.022e23 #avogadros constant, molecules/molecules
    M_to_nM = 1e9 #nano conversion, 10^9
    molecules_per_cell = (x_nanomolar*M_to_nM)/(avo_constant*cell_volume)
    return molecules_per_cell
end