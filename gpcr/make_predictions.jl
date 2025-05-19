using DifferentialEquations

function make_prediction(odeproblem, odesolverinputs, parameters)
    op = remake(odeproblem, p=parameters)
    predicted = DifferentialEquations.solve(op, odesolverinputs["solver"], abstol=odesolverinputs["abstol"], reltol=odesolverinputs["reltol"], saveat=odesolverinputs["saveat"])
    return predicted
end