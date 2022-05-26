
"""
PoissonParms

A Poisson distribution
"""
struct PoissonParms <: DistributionParmSet
  λ::Real
  useful_max::Real
end
function PoissonParms(λ)
  return PoissonParms(λ, 5λ)
end

usefulMax(d::PoissonParms) = d.useful_max

function draw(distribution::PoissonParms)
  return inv_pdf(distribution, rand())
end

function pdf(distribution::PoissonParms, t::Float64)
  return distribution.λ * exp(-distribution.λ*t)
end

function inv_pdf(distribution::PoissonParms, p::Real)
  return -log.(p/distribution.λ)/distribution.λ
end

function describe(distribution::PoissonParms)
  return "poisson λ=$(distribution.λ)"
end
