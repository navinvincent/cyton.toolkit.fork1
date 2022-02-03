using  Plots, Distributions

struct FixedDistributionParms <: DistributionParmSet
  x::Float64
  "This provides a hint for setting plot axis limits"
  useful_max::Float64
  function FixedDistributionParms(x::Float64)
    return new(x, 3x)
  end
end

function draw(distribution::FixedDistributionParms)
  return distribution.x
end

function describe(distribution::LogNormalParms)
  return "fixed at $(distribution.x)"
end

function pdf(d::FixedDistributionParms, t::Float64)
  if t == d.x
    return Inf
  else
    return 0.0
  end
end
