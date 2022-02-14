abstract type DistributionParmSet end

export DistributionParmSet, 
  FixedDistributionParms, 
  LogNormalParms, 
  PoissonParms, 
  draw,
  plotPdf

include("poisson.jl")
include("log_normal.jl")
include("fixed.jl")

function plotPdf(distribution::DistributionParmSet, max::Float64=distribution.useful_max)
  t = collect(0:0.1:max)
  p = pdf.(Ref(distribution), t)
  return Plots.plot(t, p)
end
  