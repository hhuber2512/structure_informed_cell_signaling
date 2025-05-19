using DifferentialEquations

function calculate_all_quantities(sol)
    n = length(sol.t)
    sim = sol.u

    #Defined in Yi
    total_G_protein = 10000 #molecules/cell, reported by Yi 
    active_G_protein = [sim[i][5] for i in 1:n] #Ga
    fraction_active_G_protein = active_G_protein./total_G_protein

    bound_ligand = [sim[i][3] for i in 1:n]

    return Dict("fraction_active_G" => fraction_active_G_protein, "RL" => bound_ligand)
end
