"""
Implements the Streched Cell Cycle Model:

Stretched cell cycle model for proliferating lymphocytes
Dowling, Mark R.Kan, AndreyHeinzel, SusanneZhou, Jie H. S.Marchingo, Julia M.
PROCEEDINGS OF THE NATIONAL ACADEMY OF SCIENCES 2014
https://dx.doi.org/10.1073/pnas.1322420111
"""

using Plots

import Cyton: 
  DistributionParmSet,
  draw,
  LogNormalParms,
  FateTimer,
  addTimer,
  Cell,
  step

# CpG stimulated B cells
# λ   = LogNormalParms(12.34, 3.48; natural=true)
λ   = LogNormalParms(1.2, 0.2)
kG1 = 0.27

mutable struct CycleTimer <: FateTimer
  divTime::Float64
  kG1::Float64
  startTime::Float64
end

"Constructor for new cells at t=0"
CycleTimer(λ::DistributionParmSet, kG1::Float64) = CycleTimer(draw(λ), kG1, 0)

"The step function for the cycle timer"
function step(cycle::CycleTimer, time::Float64, Δt::Float64)
  if time >= cycle.divTime + cycle.startTime
    # Cell has (fake) divided, reset the timer.
    cycle.startTime = time
  end
end

function remaining(cycle::CycleTimer, time::Float64) 
  t = cycle.divTime + cycle.startTime - time
  if t < -0.8
    println(".")
    println(cycle.divTime)
    println(cycle.startTime)
  end
  return t
end
remaining(cell::Cell, time::Float64) = remaining(cell.timers[1], time)

function stretchedCellFactory(birth::Float64=0.0)
  cell = Cell(birth)
  addTimer(cell, CycleTimer(λ, kG1))
  return cell
end

model = createModel(100, stretchedCellFactory)
model.properties[:Δt] = 0.1

function run(model::AgentBasedModel, runDuration::Float64)
  
  Δt = model.properties[:Δt]
  for _ in 0:Δt:runDuration
    step(model)
  end

  cellAgents = values(model.agents)
  Ncells = length(cellAgents)
  survivalTimes = zeros(Ncells)

  for (i, a) in enumerate(cellAgents)
    cell = a.cell
    survivalTimes[i] = remaining(cell, runDuration)
  end
  sort!(survivalTimes)
  tm = 1 .- collect(1:1:Ncells) ./ Ncells

  h = plot(survivalTimes, tm)
  display(h)

  println("Done at model time=$(model.properties[:step_cnt]*model.properties[:Δt])")
end
