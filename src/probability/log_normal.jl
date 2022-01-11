using Distributions: LogNormal

struct LogNormalParms <: DistributionParmSet
  μ::Float64
  σ::Float64
end

function draw(distribution::LogNormalParms)
  return rand(LogNormal(distribution.μ, distribution.σ))
end

function draw(distribution::LogNormalParms, n::Int64)
  return rand(LogNormal(distribution.μ, distribution.σ), n)
end

function pdf(distribution::LogNormalParms, t::Float64)
  return pdf(LogNormal(distribution.μ, distribution.σ), t)
end

function describe(distribution::LogNormalParms)
  return "log normal μ=$(distribution.μ) σ=$(distribution.σ)"
end

