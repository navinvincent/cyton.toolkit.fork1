abstract type DistributionParmSet end

export DistributionParmSet, 
  FixedDistributionParms, 
  LogNormalParms, 
  NormalParms, 
  PoissonParms, 
  draw,
  plotPdf,
  describe

draw(d::DistributionParmSet) = error("draw not implemented for $(typeof(d)).")

include("poisson.jl")
include("log_normal.jl")
include("normal.jl")
include("fixed.jl")

function plotPdf(distribution::DistributionParmSet, max::Real=distribution.useful_max)
  t = collect(0:0.1:max)
  p = pdf.(Ref(distribution), t)
  return Plots.plot(t, p)
end
  