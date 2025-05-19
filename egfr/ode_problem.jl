using DifferentialEquations

"""

"""
function egfr_ode!(du, u, p, t)
    EGF, R, Ra, R2, RP, RPLCg, RPLCgP, RG, RGS, RSh, RShP, RShG, RShGS, GS, ShP, ShG, PLCg, PLCgP, PLCgl, Grb, Shc, SOS, ShGS = u
    K16, K4, K8, V16, V4, V8, k10b, k10f, k11b, k11f, k12b, k12f, k13b, k13f, k14b, k14f, k15b, k15f, k17b, k17f, k18b, k18f, k19b, k19f, k1b, k1f, k20b, k20f, k21b, k21f, k22b, k22f, k23b, k23f, k24b, k24f, k25b, k25f, k2b, k2f, k3b, k3f, k5b, k5f, k6b, k6f, k7b, k7f, k9b, k9f = p
    v1 = k1f*R*EGF - k1b*Ra
    v2 = k2f*Ra*Ra - k2b*R2
    v3 = k3f*R2 - k3b*RP
    v4 = V4*RP/(K4+RP)
    v5 = k5f*RP*PLCg - k5b*RPLCg
    v6 = k6f*RPLCg - k6b*RPLCgP
    v7 = k7f*RPLCgP - k7b*RP*PLCgP
    v8 = V8*PLCgP/(K8+PLCgP)
    v9 = k9f*RP*Grb - k9b*RG
    v10 = k10f*RG*SOS - k10b*RGS
    v11 = k11f*RGS - k11b*RP*GS
    v12 = k12f*GS - k12b*Grb*SOS
    v13 = k13f*RP*Shc - k13b*RSh
    v14 = k14f*RSh - k14b*RShP
    v15 = k15f*RShP - k15b*ShP*RP
    v16 = V16*ShP/(K16+ShP)
    v17 = k17f*RShP*Grb - k17b*RShG
    v18 = k18f*RShG - k18b*RP*ShG
    v19 = k19f*RShG*SOS - k19b*RShGS
    v20 = k20f*RShGS - k20b*ShGS*RP
    v21 = k21f*ShP*Grb - k21b*ShG
    v22 = k22f*ShG*SOS - k22b*ShGS
    v23 = k23f*ShGS - k23b*ShP*GS
    v24 = k24f*RShP*GS - k24b*RShGS
    v25 = k25f*PLCgP - k25b*PLCgl
    du[1] = -v1
    du[2] = -v1
    du[3] = v1 - 2*v2
    du[4] = v2 + v4 - v3
    du[5] = v3 + v7 + v11 + v15 + v18 + v20 - v4 - v5 - v9 - v13
    du[6] = v5 - v6
    du[7] = v6 - v7
    du[8] = v9 - v10
    du[9] = v10 - v11
    du[10] = v13 - v14
    du[11] = v14 - v24 - v15 - v17
    du[12] = v17 - v18 - v19
    du[13] = v19 - v20 + v24
    du[14] = v11 + v23 - v12 - v24
    du[15] = v15 + v23 - v21 - v16
    du[16] = v18 + v21 - v22
    du[17] = v8 - v5
    du[18] = v7 - v8 - v25
    du[19] = v25
    du[20] = v12 - v9 - v17 - v21
    du[21] = v16 - v13
    du[22] = v12 - v10 - v19 - v22
    du[23] = v20 + v22 - v23
    nothing
end

"""

"""
function groundtruth_parameter_values()
    #from Kholodenko
    values = [340.0; 50.0; 100.0; 1.7; 450.0; 1.0; 0.06; 0.01; 0.0045; 0.03; 0.0001; 0.0015; 0.6; 0.09; 0.06; 6.0; 0.0009; 0.3; 0.1; 0.003; 0.0009; 0.3; 0.0214; 0.01; 0.06; 0.003; 0.00024; 0.12; 0.1; 0.003; 0.064; 0.03; 0.021; 0.1; 0.0429; 0.009; 0.03; 1.0; 0.1; 0.01; 0.01; 1.0; 0.2; 0.06; 0.05; 1.0; 0.006; 0.3; 0.05; 0.003]
    return values 
end

"""

"""
function initial_conditions(egf_dose)
    #from Table II of Kholodenko
    if egf_dose == 20.0
        u0 = [680.0; 100.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 105.0; 0.0; 0.0; 85.0; 150.0; 34.0; 0.0]
    elseif egf_dose == 2.0
        u0 = [68.0; 100.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 105.0; 0.0; 0.0; 85.0; 150.0; 34.0; 0.0]
    elseif egf_dose == 0.2
        u0 = [6.8; 100.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 105.0; 0.0; 0.0; 85.0; 150.0; 34.0; 0.0]
    else 
        print("check egf dose")
    end
    return u0
end

"""

return\\_ode\\_problem\\_default\\_inputs() \n

Function returns four inputs needed to define an ODEProblem from the DifferentialEquations package for EGFR signaling. \n
Species units are molecules, rate constants are 1/sec or 1/sec*nM \n
Default initial conditions and parameter values were reported in Kholodenko et al. \n

Should return: \n
odesys: DifferentialEquations.ODESystem \n
u0: Vector{Pair{Num, Float64}} \n
tspan: Tuple{Int64, Int64} \n
p: Vector{Pair{Num, Float64}} \n

"""

function return_ode_problem_default_inputs(egf_dose)
    odesys = egfr_ode!
    tspan = (0,120)
    p = groundtruth_parameter_values()
    u0 = initial_conditions(egf_dose)
    return odesys, u0, tspan, p
end

"""

return\\_ode\\_problem\\_solver\\_default\\_inputs()

Function returns four inputs needed to solve an ODEProblem from the DifferentialEquations package for EGFR signaling. \n
Species units are molecules, rate constants are 1/sec or 1/sec*nM \n
Default initial conditions and parameter values were reported in Kholodenko et al. \n

Should return: Dict with following keys: \n
solver: ModelingToolkit.ODESystem \n
abstol: Vector{Pair{Num, Float64}} \n
reltol: Tuple{Int64, Int64} \n
saveat: Vector{Pair{Num, Float64}} \n

"""
function return_ode_problem_solver_default_inputs()
    #absolute solver tolerance based on protein concentration deemed insignificant, 1e-5 nM, or less than 1 protein per cell
    abstol=1e-5
    #relative solver tolerance based on number of significant digits for protein concentration (5), plus 1: 1e-(5+1) 
    reltol=1e-6
    #stiff solver
    solver = QNDF()
    #all species observed at the same timepoints, can take any save_at
    saveat = deserialize("outputs/000_processed_grb_egfr_20.dict")["save_at"]
    return Dict("solver"=>solver, "abstol" => abstol, "reltol" => reltol, "saveat" => saveat)
end

function return_finegrain_saveat()
    return collect(range(start=0.0,stop=120.0, step=1))
end
