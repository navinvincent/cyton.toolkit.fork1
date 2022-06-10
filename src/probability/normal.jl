using  Plots
import Distributions as dst

"""
NormalParms

A normal distribution.
"""
struct NormalParms <: DistributionParmSet
  μ::Real
  σ::Real
  "This provides a hint for setting plot axis limits"
  useful_max::Real
end
NormalParms(μ::Real, σ::Real) = NormalParms(μ, σ, μ + 3σ)

usefulMax(d::NormalParms) = d.useful_max

sample(distribution::NormalParms) = dst.rand(dst.Normal(distribution.μ, distribution.σ))

sample(distribution::NormalParms, n::Int64) = dst.rand(dst.Normal(distribution.μ, distribution.σ), n)

pdf(distribution::NormalParms, t::Float64) = dst.pdf(dst.Normal(distribution.μ, distribution.σ), t)

describe(distribution::NormalParms) = "normal μ=$(distribution.μ) σ=$(distribution.σ)"

