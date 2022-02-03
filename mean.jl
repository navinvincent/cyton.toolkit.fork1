using Distributions, Plots

# Parameters from the Cyton2 paper, fig 5A1
λ_firstDivision = LogNormal(log(39.89), 0.28)
λ_subsequentDivision = 9.21
λ_divisionDestiny = LogNormal(log(1171.86), 0.11)
λ_lifetime = LogNormal(log(11116.8), 0.85)

Δt = 0.5
tm = 0:Δt:150
N0 = 1000

divisionSurvival = ccdf(λ_firstDivision, tm)
deathSurvival    = ccdf(λ_lifetime, tm)
destinySurvival  = ccdf(λ_divisionDestiny, tm)
destinyPdf       = pdf.(λ_divisionDestiny, tm)

# Gen 0
p_div = divisionSurvival .* destinySurvival
p_dest = cumsum([p * s for (p, s) in zip(destinyPdf, divisionSurvival)]) .* Δt
nCells0 = N0 .* deathSurvival .* (p_dest .+ p_div)

# Gen 1
gen = 1
p_div_diff = cdf(λ_firstDivision, tm .- (gen-1)*λ_subsequentDivision) - cdf(λ_firstDivision, tm .- gen*λ_subsequentDivision)
p_div = destinySurvival .* p_div_diff
p_dest = cumsum([p * s for (p, s) in zip(destinyPdf, p_div_diff)]) .* Δt
nCells1 = 2 .* N0 .* deathSurvival .* (p_dest .+ p_div)

plot(tm, [nCells0, nCells1])

