include("michelle.jl")
include("plotting.jl")

# @threads 
for parm in [parameters[1]]
  r = runModel(parm)
  bigPlot(parm, r)
end
