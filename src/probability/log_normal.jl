using Distributions, Plots

struct LogNormalParms <: DistributionParmSet
  μ::Float64
  σ::Float64
  useful_max::Float64
  function LogNormalParms(μ::Float64, σ::Float64)
    return new(μ, σ, μ+10σ)
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

function plot_pdf(distribution::LogNormalParms, max::Float64=distribution.useful_max)
  t = collect(0:0.1:max)
  p = pdf.(Ref(distribution), t)
  return Plots.plot(t, p)
end
