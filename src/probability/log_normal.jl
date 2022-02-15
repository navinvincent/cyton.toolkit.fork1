using  Plots, Distributions

struct LogNormalParms <: DistributionParmSet
  μ::Real
  σ::Real
  "This provides a hint for setting plot axis limits"
  useful_max::Real
end
function LogNormalParms(μ::Real, σ::Real; natural=false)
  if natural
    f = 1+ (σ/μ)^2
    μ_L = log(μ/sqrt(f))
    σ_L = sqrt(log(f))
    return LogNormalParms(μ_L, σ_L, exp(μ_L + 3σ_L))
  else
    return LogNormalParms(μ, σ, exp(μ + 3σ))
  end
end

function draw(distribution::LogNormalParms)
  return rand(LogNormal(distribution.μ, distribution.σ))
end

function draw(distribution::LogNormalParms, n::Int64)
  return rand(LogNormal(distribution.μ, distribution.σ), n)
end

function pdf(distribution::LogNormalParms, t::Float64)
  return Distributions.pdf(LogNormal(distribution.μ, distribution.σ), t)
end

function describe(distribution::LogNormalParms)
  return "log normal μ=$(distribution.μ) σ=$(distribution.σ)"
end
