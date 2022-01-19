using  Plots, Distributions

struct LogNormalParms <: DistributionParmSet
  μ::Float64
  σ::Float64
  "This provides a hint for setting plot axis limits"
  useful_max::Float64
  function LogNormalParms(μ::Float64, σ::Float64)
    return new(μ, σ, μ+10σ)
  end
  function LogNormalParms(μ::Float64, σ::Float64, useful_max::Float64)
    return new(μ, σ, useful_max)
  end
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
