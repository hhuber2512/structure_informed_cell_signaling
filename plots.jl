
function egfr_return_kon_indices()
    return [10,14,16,22,26,28,32,46,52].-2
end

function egfr_return_koff_indices()
    return [9,13,15,21,25,27,31,45,51].-2
end

function egfr_return_nonbinding_indices()
    n_parameters = 50
    binding_indices = cat(egfr_return_koff_indices(), egfr_return_kon_indices(), dims=1)
    n_binding_parameters = length(binding_indices)
    return collect(1:n_parameters)[[all(fill(collect(1:n_parameters)[i],n_binding_parameters) .!= binding_indices) for i in 1:n_parameters]]
end

function gpcr_return_kon_indices()
    return [1]
end

function gpcr_return_koff_indices()
    return [2]
end

function gpcr_return_nonbinding_indices()
    return [3,4,5,6,7,8]
end

function return_gpcr_parameters()
    return ["k1f", "k1b", "k2", "k3", "k4", "k5", "k6", "k7"]
end

function convert_inverse_nM_molecules(x_nanomolar)
    #N.B. converting from 1/(nM*sec) to 1/(molecules per cell * sec)
    cell_volume = 8.6e-14 #Liters, average cell volume from bionumbers - https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=4&id=111978
    avo_constant = 6.022e23 #avogadros constant, molecules/molecules
    M_to_nM = 1e9 #nano conversion, 10^9
    molecules_per_cell = (x_nanomolar*M_to_nM)/(avo_constant*cell_volume)
    return molecules_per_cell
end

function return_gpcr_parameter_values()
    return [convert_inverse_nM_molecules(2.0e6/1e9); 0.01; 1.0; 1.0E-5; 4.0; 4.0E-4; 0.0040; 0.11]
end

function return_egfr_parameters()
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

function return_egfr_parameter_values()
    values = [340.0; 50.0; 100.0; 1.7; 450.0; 1.0; 0.06; 0.01; 0.0045; 0.03; 0.0001; 0.0015; 0.6; 0.09; 0.06; 6.0; 0.0009; 0.3; 0.1; 0.003; 0.0009; 0.3; 0.0214; 0.01; 0.06; 0.003; 0.00024; 0.12; 0.1; 0.003; 0.064; 0.03; 0.021; 0.1; 0.0429; 0.009; 0.03; 1.0; 0.1; 0.01; 0.01; 1.0; 0.2; 0.06; 0.05; 1.0; 0.006; 0.3; 0.05; 0.003]
    return values
end

#plotting hyperparameters
function return_plot_inputs(case)
    if case == "p_shc_20"
        x = deserialize("egfr/outputs/000_processed_$(case).dict")["save_at"]
        response = deserialize("egfr/outputs/000_processed_$(case).dict")["response"]
        avg_error = deserialize("egfr/outputs/000_processed_$(case).dict")["average_error"]
        x_finegrain = collect(range(start=0.0,stop=120.0, step=1))
        ylabel = "Percent phosphorylated SHC, 20nM EGF"
        xlabel = "Time (sec)"
        xlims = [0,125]
        ylims = [0,105]
    elseif case == "p_egfr_20"
        x = deserialize("egfr/outputs/000_processed_$(case).dict")["save_at"]
        response = deserialize("egfr/outputs/000_processed_$(case).dict")["response"]
        avg_error = deserialize("egfr/outputs/000_processed_$(case).dict")["average_error"]
        x_finegrain = collect(range(start=0.0,stop=120.0, step=1))
        ylabel = "Percent phosphorylated EGFR, 20nM EGF"
        xlabel = "Time (sec)"
        xlims = [0,125]
        ylims = [0,105]
    elseif case == "p_egfr_2"
        x = deserialize("egfr/outputs/000_processed_$(case).dict")["save_at"]
        response = deserialize("egfr/outputs/000_processed_$(case).dict")["response"]
        avg_error = deserialize("egfr/outputs/000_processed_$(case).dict")["average_error"]
        x_finegrain = collect(range(start=0.0,stop=120.0, step=1))
        ylabel = "Percent phosphorylated EGFR, 2nM EGF"
        xlabel = "Time (sec)"
        xlims = [0,125]
        ylims = [0,105]
    elseif case == "p_egfr_02"
        x = deserialize("egfr/outputs/000_processed_$(case).dict")["save_at"]
        response = deserialize("egfr/outputs/000_processed_$(case).dict")["response"]
        avg_error = deserialize("egfr/outputs/000_processed_$(case).dict")["average_error"]
        x_finegrain = collect(range(start=0.0,stop=120.0, step=1))
        ylabel = "Percent phosphorylated EGFR, 0.2nM EGF"
        xlabel = "Time (sec)"
        xlims = [0,125]
        ylims = [0,105]
    elseif case == "p_plcg_20"
        x = deserialize("egfr/outputs/000_processed_$(case).dict")["save_at"]
        response = deserialize("egfr/outputs/000_processed_$(case).dict")["response"]
        avg_error = deserialize("egfr/outputs/000_processed_$(case).dict")["average_error"]
        x_finegrain = collect(range(start=0.0,stop=120.0, step=1))
        ylabel = "Percent phosphorylated PLCg, 20 nM EGF"
        xlabel = "Time (sec)"
        xlims = [0,125]
        ylims = [0,15]
    elseif case == "p_plcg_2"
        x = deserialize("egfr/outputs/000_processed_$(case).dict")["save_at"]
        response = deserialize("egfr/outputs/000_processed_$(case).dict")["response"]
        avg_error = deserialize("egfr/outputs/000_processed_$(case).dict")["average_error"]
        x_finegrain = collect(range(start=0.0,stop=120.0, step=1))
        ylabel = "Percent phosphorylated PLCg, 2nM EGF"
        xlabel = "Time (sec)"
        xlims = [0,125]
        ylims = [0,15]
    elseif case == "grb_shc_20"
        x = deserialize("egfr/outputs/000_processed_$(case).dict")["save_at"]
        response = deserialize("egfr/outputs/000_processed_$(case).dict")["response"]
        avg_error = deserialize("egfr/outputs/000_processed_$(case).dict")["average_error"]
        x_finegrain = collect(range(start=0.0,stop=120.0, step=1))
        ylabel = "Percent GRB2 in SHC IP, 20nM EGF"
        xlabel = "Time (sec)"
        xlims = [0,125]
        ylims = [0,45]
    elseif case == "grb_egfr_20"
        x = deserialize("egfr/outputs/000_processed_$(case).dict")["save_at"]
        response = deserialize("egfr/outputs/000_processed_$(case).dict")["response"]
        avg_error = deserialize("egfr/outputs/000_processed_$(case).dict")["average_error"]
        x_finegrain = collect(range(start=0.0,stop=120.0, step=1))
        ylabel = "Percent GRB2 in EGFR IP, 20nM EGF"
        xlabel = "Time (sec)"
        xlims = [0,125]
        ylims = [0,25]
    elseif case == "active_g_dose_response"
        x = log10.(deserialize("gpcr/outputs/000_processed_$(case).dict")["ligand_stimulation (nM)"])
        response = deserialize("gpcr/outputs/000_processed_$(case).dict")["response"]
        avg_error = deserialize("gpcr/outputs/000_processed_$(case).dict")["average_error"]
        x_finegrain = collect(range(start=-2,stop=3,length=50))
        ylabel = "G-Protein Activation Response"
        xlabel = "log[Alpha-Factor](nM)"
        xlims = [-2,3]
        ylims = [0,1.2]
    elseif case == "rl_dose_response"
        x = log10.(deserialize("gpcr/outputs/000_processed_$(case).dict")["ligand_stimulation (nM)"])
        response = deserialize("gpcr/outputs/000_processed_$(case).dict")["response"]
        avg_error = deserialize("gpcr/outputs/000_processed_$(case).dict")["average_error"]
        x_finegrain = collect(range(start=-2,stop=3,length=50))
        ylabel = "RL Response"
        xlabel = "log[Alpha-Factor](nM)"
        xlims = [-2,3]
        ylims = [0,1.2]
    elseif case == "active_g_timecourse"
        x = deserialize("gpcr/outputs/000_processed_$(case).dict")["save_at"]
        response = deserialize("gpcr/outputs/000_processed_$(case).dict")["response"]
        avg_error = deserialize("gpcr/outputs/000_processed_$(case).dict")["average_error"]
        x_finegrain = collect(range(start=0.0, stop=600.0, length=100))
        ylabel = "Fraction Active G Protein"
        xlabel = "Time (sec)"
        xlims = [0,650]
        ylims = [0,0.5]
    end
    return ylims, xlims, xlabel, ylabel, x_finegrain,avg_error, response, x
end

function plot_quantiles(ylims, xlims, xlabel, ylabel, x_finegrain, avg_error, response, x, color, quantiles_reg, quantiles_unreg)
    marker_size = 7
    linewidth = 3
    p = plot(x_finegrain, quantiles_unreg["lower_quantile"], linewidth=linewidth; fillrange=quantiles_unreg["upper_quantile"], label="baseline", color=color, fillalpha=0.5)
    plot!(x_finegrain, quantiles_unreg["upper_quantile"], label=false, color=color, linewidth=linewidth)
    plot!(x_finegrain, quantiles_unreg["median"], label="median baseline", color=color, linestyle = :dot, linewidth=linewidth)
    p = plot!(x_finegrain, quantiles_reg["lower_quantile"], linewidth=linewidth; fillrange=quantiles_reg["upper_quantile"], color=color, fillstyle = :|,label="augmented")
    plot!(x_finegrain, quantiles_reg["upper_quantile"], label=false, color=color, linewidth=linewidth)
    plot!(x_finegrain, quantiles_reg["median"], label="median augmented", color=color, linestyle = :dash,linewidth=linewidth)
    scatter!(x, response, label="experimental data", color="black", yerr=avg_error, markersize=marker_size)
    xlabel!(xlabel)
    ylabel!(ylabel)
    xlims!(xlims[1],xlims[2])
    ylims!(ylims[1],ylims[2])
    plot!(dpi=300,size=(690,690),legendfontsize=8,titlefontsize=8,
    tickfontsize=8,guidefontsize=8,left_margin = 8Plots.mm, bottom_margin = 5Plots.mm)
    return p
end

function plot_quantiles_no_experimental_data(my_species, x_finegrain, xlims, color, quantiles_reg, quantiles_unreg)
    p = plot(x_finegrain, quantiles_unreg["lower_quantile"]; fillrange=quantiles_unreg["upper_quantile"], label="baseline", color=color, fillalpha=0.5)
    plot!(x_finegrain, quantiles_unreg["upper_quantile"], label=false, color=color)
    plot!(x_finegrain, quantiles_unreg["median"], label="median baseline", color=color, linestyle = :dot, linewidth=2)
    p = plot!(x_finegrain, quantiles_reg["lower_quantile"]; fillrange=quantiles_reg["upper_quantile"], color=color, fillstyle = :|,label="augmented")
    plot!(x_finegrain, quantiles_reg["upper_quantile"], label=false, color=color)
    plot!(x_finegrain, quantiles_reg["median"], label="median augmented", color=color, linestyle = :dash,linewidth=2)
    ylabel!(my_species)
    xlabel!("Time (sec)")
    xlims!(xlims[1], xlims[2])
    plot!(dpi=300,size=(690,690),legendfontsize=8,titlefontsize=8,
    tickfontsize=8,guidefontsize=8,left_margin = 8Plots.mm, bottom_margin = 5Plots.mm)
    return p
end

function return_egfr_species()
    species_ordered = ["EGF", "R", "Ra", "R2", "RP", "RPLCg", "RPLCgP", "RG", "RGS", "RSh", "RShP", "RShG", "RShGS", "GS", "ShP","ShG","PLCg","PLCgP", "PLCgl","Grb","Shc","SOS","ShGS"]
    return species_ordered
end

function return_gpcr_species()
    species_ordered = ["R", "L", "RL", "G", "Ga"]
    return species_ordered
end

function return_egfr_ligand_doses_naming()
    return ["20", "2", "02"]
end

function return_gpcr_ligand_doses_naming()
    return ["01", "1", "2", "5", "10", "20", "50", "100", "1000"]
end

function sensitivity_scatter_plot(sensitivity, diff_med, color)
    p = scatter(sensitivity, diff_med, color=color,clims=(-1,1))
    plot!(dpi=300,size=(690,690),legendfontsize=8,titlefontsize=8,
    tickfontsize=8,guidefontsize=8,left_margin = 8Plots.mm, bottom_margin = 5Plots.mm)
    return p
end