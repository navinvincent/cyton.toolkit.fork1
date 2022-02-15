"""
Implements the Streched Cell Cycle Model:

Stretched cell cycle model for proliferating lymphocytes
Dowling, Mark R.Kan, AndreyHeinzel, SusanneZhou, Jie H. S.Marchingo, Julia M.
PROCEEDINGS OF THE NATIONAL ACADEMY OF SCIENCES 2014
https://dx.doi.org/10.1073/pnas.1322420111
"""

using Agents
import Gadfly
# plotlyjs()
import Base: step # work around bug in adding additional methods with this name

using Cyton

@enum Phase G1 S G2M
# CpG stimulated B cells
λ   = LogNormalParms(12.34, 3.48; natural=true)
kG1 = 0.27
kS  = 0.57
# BRDU labelled and unlabelled floursence level
brduLo = LogNormalParms(log(250), log(2))
brduHi = LogNormalParms(log(7000), log(2))
# 7AAD DNA 2x and 4x labelling
dnaLo = NormalParms(75000, 5000)
dnaHi = NormalParms(150000, 10000)

mutable struct CycleTimer <: FateTimer
  divTime::Float64
  kG1::Float64
  kS::Float64
  startTime::Float64
  brdu::Float64
end

"Constructor for new cells at t=0"
function CycleTimer(λ::DistributionParmSet, kG1::Float64)
  # Timers are desynchronised by distributing them randomly 
  # over their cycle.
  l = draw(λ)
  s = - l * rand()
  CycleTimer(l, kG1, kS, s, draw(brduLo))
end

function phase(cycle::CycleTimer, time::Float64)
  δt = time - cycle.startTime
  divTime = cycle.divTime
  kG1 = cycle.kG1 * divTime
  kS = cycle.kS * divTime
  if δt <= kG1
    return G1
  end
  if δt <= kG1 + kS
    return S
  end
  return G2M
end

"The step function for the cycle timer"
function step(cycle::CycleTimer, time::Float64, Δt::Float64)
  if time >= cycle.divTime + cycle.startTime
    # Cell has (fake) divided, reset the timer.
    cycle.startTime = time
  end
end

remaining(cycle::CycleTimer, time::Float64) = cycle.divTime + cycle.startTime - time
remaining(cell::Cell, time::Float64) = remaining(cell.timers[1], time)

function stretchedCellFactory(birth::Float64=0.0)
  cell = Cell(birth)
  addTimer(cell, CycleTimer(λ, kG1))
  return cell
end

function runModel!(model::AgentBasedModel, runDuration::Float64, stimulus::Stimulus, callback::Function=(_) -> nothing)
  Δt = model.properties[:Δt]
  for _ in 0:Δt:runDuration
    step(model, stimulus)
    callback(model)
  end
end

function survivalCurves(model::AgentBasedModel)
  cellAgents = values(model.agents)
  totalAlpha = Vector{Float64}()
  g1Alpha = Vector{Float64}()
  sg2mAlpha = Vector{Float64}()
  
  runDuration = model_time(model)
  for a in cellAgents
    cell = a.cell
    push!(totalAlpha, remaining(cell, runDuration))
    cycle = cell.timers[1]
    if phase(cycle, runDuration) == G1
      r = cycle.divTime * cycle.kG1 + cycle.startTime - runDuration
      push!(g1Alpha, r)
    else
      r = (cycle.divTime + cycle.startTime - runDuration) * (1 - cycle.kG1)
      push!(sg2mAlpha, r)
    end
  end

  sort!(totalAlpha)
  n = length(totalAlpha)
  tmT = 1 .- collect(1:1:n) ./ n
  sort!(g1Alpha)
  n = length(g1Alpha)
  tmG1 = 1 .- collect(1:n)/n
  sort!(sg2mAlpha)
  n = length(sg2mAlpha)
  tmSg2m = 1 .- collect(1:n)/n

  h = plot()
  plot!(g1Alpha, tmG1, label="G1", lw=3, lc="green")
  plot!(sg2mAlpha, tmSg2m, label="S/G2/M", lw=3, lc="red")
  plot!(totalAlpha, tmT, label="Total", lw=3, lc="blue")
  xlabel!("Time (h)")
  ylabel!("proportion")
  yaxis!(:log)
  yticks!([10.0^x for x in [-3, -2, -1, 0]])
  display(h)

  println("Done at model time: $(model.properties[:step_cnt]*model.properties[:Δt])")
end

function timeInS(cycle::CycleTimer, time::Float64)
  δt = time - cycle.startTime
  divTime = cycle.divTime
  kG1 = cycle.kG1 * divTime
  kS = cycle.kS * divTime
  if δt <= kG1
    error("Cell is not in S")
  end
  if δt <= kG1 + kS
    return (δt-kG1)/kS
  end
  error("Cell is not in S")
end

function dnaStainLevels(model::AgentBasedModel)
  time = model_time(model)
  NCells = length(model.agents)
  brdus = zeros(NCells)
  dnas = zeros(NCells)
  for (i, a) in enumerate(values(model.agents))
    c = a.cell
    fate = c.timers[1]
    brdus[i] = fate.brdu
    p = phase(fate, time)
    
    l = draw(dnaLo)
    if p == G1
      dnas[i] = l
      continue
    end

    h = draw(dnaHi)
    if p == S
      t = timeInS(fate, time)
      dnas[i] = l + (h-l)*t
      continue
    end
    
    if p == G2M
      dnas[i] = h
    end
  end
  return (dnas, brdus)
end

struct BrduStimulus <: Stimulus
  pulseStart::Float64
  pulseEnd::Float64
end

function Cyton.stimulate(cell::Cell, stim::BrduStimulus, time::Float64)
  cycle = cell.timers[1]
  if (stim.pulseStart ≤ time ≤ stim.pulseEnd) && (phase(cycle, time) == S)
    cycle.brdu = draw(brduHi)
  end
end
stim = BrduStimulus(0, 0.5)

println("-------------------- start --------------------")
model = createModel(20000, stretchedCellFactory)
model.properties[:Δt] = 0.01

rt = @timed runModel!(model, 1.0, stim)
println("Elaspsed time: $(rt[2])")

(dnas, brdus) = dnaStainLevels(model);;;;;;;;;
# survivalCurves(model)

plot(x=dnas, 
y=brdus, 
Geom.histogram2d, 
Guide.xlabel("DNA"), 
Guide.ylabel("BRDU"), 
Coord.cartesian(xmin=0, xmax=200000, ymin=1, ymax=5), 
Scale.y_log10, 
Theme(background_color="white"))

# p |> PNG("sin_pi.png", 6inch, 4inch)