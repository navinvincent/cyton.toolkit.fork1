module DistributionParms

export DistributionParmSet, PoissonParms, LogNormalParms, draw, plot_pdf

abstract type DistributionParmSet end

include("poisson.jl")
include("log_normal.jl")

end
