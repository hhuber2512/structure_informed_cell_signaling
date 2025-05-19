

"""

"""
function return_index_order_of_data_for_likelihood_w_pshc_test(egf_dosage)
    if egf_dosage == 20.0
        return ["p_egfr", "p_plcg", "grb_egfr", "grb_shc"]
    elseif egf_dosage == 2.0
        return ["p_egfr", "p_plcg"]
    elseif egf_dosage == 0.2
        return ["p_egfr"]
    end
end

"""

"""
function return_index_order_of_data_for_likelihood_w_pshc_grb2shc_test(egf_dosage)
    if egf_dosage == 20.0
        return ["p_egfr", "p_plcg", "grb_egfr"]
    elseif egf_dosage == 2.0
        return ["p_egfr", "p_plcg"]
    elseif egf_dosage == 0.2
        return ["p_egfr"]
    end
end

"""

"""
function return_index_order_of_data_for_likelihood_w_egfr_test(egf_dosage)
    if egf_dosage == 20.0
        return ["p_shc", "p_plcg", "grb_egfr", "grb_shc"]
    elseif egf_dosage == 2.0
        return ["p_plcg"]
    end
end

"""

"""
function likelihood_w_pshc_test(experimental_quantities, egf_dosage)
    quantities_per_dose = Array{Float64}(undef,0)
    species_index = return_index_order_of_data_for_likelihood_w_pshc_test(egf_dosage)
    [append!(quantities_per_dose, experimental_quantities[j]) for j in species_index]
    return quantities_per_dose
end

"""

"""
function likelihood_w_egfr_test(experimental_quantities, egf_dosage)
    quantities_per_dose = Array{Float64}(undef,0)
    species_index = return_index_order_of_data_for_likelihood_w_egfr_test(egf_dosage)
    [append!(quantities_per_dose, experimental_quantities[j]) for j in species_index]
    return quantities_per_dose
end

"""

"""
function likelihood_w_pshc_grb2shc_test(experimental_quantities, egf_dosage)
    quantities_per_dose = Array{Float64}(undef,0)
    species_index = return_index_order_of_data_for_likelihood_w_pshc_grb2shc_test(egf_dosage)
    [append!(quantities_per_dose, experimental_quantities[j]) for j in species_index]
    return quantities_per_dose
end

"""

"""
function return_ligand_dose_order_for_likelihood_w_pshc_test()
    return [20.0, 2.0, 0.2]
end

"""

"""
function return_ligand_dose_order_for_likelihood_w_pshc_grb2shc_test()
    return [20.0, 2.0, 0.2]
end

"""

"""
function return_ligand_dose_order_for_likelihood_w_egfr_test()
    return [20.0, 2.0]
end