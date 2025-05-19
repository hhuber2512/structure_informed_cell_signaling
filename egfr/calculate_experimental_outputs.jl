using DifferentialEquations
"""

"""
function calculate_all_quantities(sol)
    n = length(sol.t)
    sim = sol.u

    #defined in Kholodenko
    total_phosphorylated_EGFR = [sum(sim[i][[9, 5, 7, 8, 11, 13, 12, 10, 6]]) for i in 1:n]
    total_EGFR = [sum(sim[i][[3, 9, 5, 7, 8, 11, 13, 12, 10, 6, 2, 4]]) for i in 1:n]
    percent_phosphorylated_EGFR = total_phosphorylated_EGFR./total_EGFR
    
    total_phosphorylated_SHC = [sum(sim[i][[23, 11, 13, 12, 15, 16]]) for i in 1:n]
    total_SHC = [sum(sim[i][[21, 23, 11, 13, 12, 15, 10, 16]]) for i in 1:n]
    percent_phosphorylated_SHC = total_phosphorylated_SHC./total_SHC

    total_phosphorylated_PLCg = [sum(sim[i][[7,18]]) for i in 1:n]
    total_PLCg = [sum(sim[i][[ 7, 17, 19, 6, 18]]) for i in 1:n]
    percent_phosphorylated_PLCg = total_phosphorylated_PLCg./total_PLCg
    
    GRB2_EGFR = [sum(sim[i][[ 9, 8, 13, 12]]) for i in 1:n]
    percent_GRB2_EGFR = GRB2_EGFR./total_EGFR

    GRB2_SHC = [sum(sim[i][[23, 13, 12, 16]]) for i in 1:n]
    percent_GRB2_SHC = GRB2_SHC./total_SHC

    return Dict("p_egfr" => 100*percent_phosphorylated_EGFR, "p_shc" => 100*percent_phosphorylated_SHC, 
    "p_plcg" => 100*percent_phosphorylated_PLCg, "grb_egfr" => 100*percent_GRB2_EGFR, "grb_shc" => 100*percent_GRB2_SHC)
end


function return_ligand_dose_order_for_naming()
    return ["20", "2", "02"]
end