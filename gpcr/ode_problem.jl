using DifferentialEquations
include("misc.jl")

function gpcr_ode!(du, u, p, t)
    Gt = 10000 #molecules, given by Yi et al
    R, L, RL, G, Ga = u
    Gd = Gt - G - Ga
    Gbg = Gt - G
    k1f, k1b, k2, k3, k4, k5, k6, k7 = p
    v1 = L*R*k1f
    v2 = RL*k1b
    v3 = Gd*Gbg*k2
    v4 = RL*G*k3
    v5 = k4
    v6 = R*k5
    v7 = RL*k6
    v8 = Ga*k7
    du[1] = -v1 + v2 + v5 - v6
    du[2] = 0 #constant in Yi el al.
    du[3] = v1 - v2 - v7
    du[4] = v3 - v4 
    du[5] = v4 - v8
    nothing
end

function groundtruth_parameter_values()
    #from Yi et al.
    #forward binding rate of receptor/ligand given in 1/(M*sec)
    #all other units are molecules per cell and seconds
    #thus, must convert ligand concentration in M to molecules per cell
    #note, function converts from 1/(nM*sec), so we convert to 1/nM from 1/M first
    values = [convert_inverse_nM_molecules(2.0e6/1e9); 0.01; 1.0; 1.0E-5; 4.0; 4.0E-4; 0.0040; 0.11]
    return values 
end

function return_ligand_dose_nM()
    return [0.1, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
end

function normalize_to_dose()
    return 1000.0
end

"""

"""
function initial_conditions(ligand_dose)
    #from Yi et al.
    #ligand dose given in nM, per Yi et al.
    #all other units are molecules per cell
    #thus, must convert ligand concentration in nM to molecules per cell
    u0 = [10000.0; convert_nM_molecules(ligand_dose); 0.0; 10000.0; 0.0]
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

function return_ode_problem_default_inputs(ligand_dose)
    odesys = gpcr_ode!
    tspan = (0,600)
    p = groundtruth_parameter_values()
    u0 = initial_conditions(ligand_dose)
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
function return_ode_problem_solver_default_inputs(output)
    #absolute solver tolerance based on protein concentration deemed insignificant, 1e-5 nM, or less than 1 protein per cell
    #made a bit smaller, due to instability at low ligand concentrations
    abstol=1e-7
    #relative solver tolerance based on number of significant digits for protein concentration (5), plus 1: 1e-(5+1) 
    #made a bit smaller, due to instability at low ligand concentrations
    reltol=1e-8
    #stiff solver
    solver = QNDF()
    #save at changes based on whether we're looking at dose response or timecourse output
    #times taken from Yi et al.
    if output == "timecourse"
        saveat = [0, 10, 30, 60, 120, 210, 300, 450, 600]
    elseif output == "dose response"
        saveat = [60] #needs to be vector
    end
    return Dict("solver"=>solver, "abstol" => abstol, "reltol" => reltol, "saveat" => saveat)
end

function return_finegrain_saveat()
    return collect(range(start=0.0, stop=600.0, length=100))
end

function return_finegrain_ligand_dose_nM()
    #note, for plotting, returns on nM scale - must be converted for simulation
    #note, range taken on log10(nM) scale, to match x axis of plots
    #thus, need to take range as exponent of 10
    return 10.0.^collect(range(start=-2,stop=3, length=50))
end
