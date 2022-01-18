
struct PoissonParms <: DistributionParmSet
  λ::Float64
  useful_max::Float64
  function PoissonParms(λ)
    return new(λ, 5λ)
  end
end

function draw(distribution::PoissonParms)
  return inv_pdf(distribution, rand())
end

function pdf(distribution::PoissonParms, t::Float64)
  return distribution.λ * exp(-distribution.λ*t)
end

function inv_pdf(distribution::PoissonParms, p::Float64)
  return -log.(p/distribution.λ)/distribution.λ
end

function describe(distribution::PoissonParms)
  return "poisson λ=$(distribution.λ)"
end
