using  Plots
import Distributions as dst

struct NormalParms <: DistributionParmSet
  μ::Real
  σ::Real
  "This provides a hint for setting plot axis limits"
  useful_max::Real
end
NormalParms(μ::Real, σ::Real) = NormalParms(μ, σ, μ + 3σ)


function draw(distribution::NormalParms)
  return dst.rand(dst.Normal(distribution.μ, distribution.σ))
end

function draw(distribution::NormalParms, n::Int64)
  return dst.rand(dst.Normal(distribution.μ, distribution.σ), n)
end

function pdf(distribution::NormalParms, t::Float64)
  return dst.pdf(dst.Normal(distribution.μ, distribution.σ), t)
end

function describe(distribution::NormalParms)
  return "normal μ=$(distribution.μ) σ=$(distribution.σ)"
end
