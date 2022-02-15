using  Plots, Distributions

struct FixedDistributionParms <: DistributionParmSet
  x::Real
  "This provides a hint for setting plot axis limits"
  useful_max::Real
  function FixedDistributionParms(x::Real)
    return new(x, 3x)
  end
end

function draw(distribution::FixedDistributionParms)
  return distribution.x
end

function describe(distribution::LogNormalParms)
  return "fixed at $(distribution.x)"
end

function pdf(d::FixedDistributionParms, t::Real)
  if t == d.x
    return Inf
  else
    return 0.0
  end
end
