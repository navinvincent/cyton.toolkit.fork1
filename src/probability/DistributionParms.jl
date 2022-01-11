module DistributionParms

export DistributionParmSet, PoissonParms, LogNormalParms, draw

abstract type DistributionParmSet end

include("poisson.jl")
include("log_normal.jl")

end
