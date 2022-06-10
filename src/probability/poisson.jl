
"""
PoissonParms

A Poisson distribution
"""
struct PoissonParms <: DistributionParmSet
  λ::Real
  useful_max::Real
end
PoissonParms(λ) = PoissonParms(λ, 5λ)

usefulMax(d::PoissonParms) = d.useful_max

sample(distribution::PoissonParms) = inv_pdf(distribution, rand())

pdf(distribution::PoissonParms, t::Float64) = distribution.λ * exp(-distribution.λ*t)

inv_pdf(distribution::PoissonParms, p::Real) = -log.(p/distribution.λ)/distribution.λ

describe(distribution::PoissonParms) = "poisson λ=$(distribution.λ)"

