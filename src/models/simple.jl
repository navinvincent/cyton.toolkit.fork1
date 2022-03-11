
"""
A simple cell has the following behaviour
- It does not divide prior to stimulation
- After stimulation there is a delay to first division
- After first division the cell will divide until it reaches destiny
- Most cells will die after a period of time.
- A proportion of cells will not die - memory cells
  - memory cell fate occurs at constant (slow) rate
  - memory cells don't divide

In this simple, proof of concept model:
- all lifetimes are drawn from LogNormalParms
- Time to next division are initiated de novo on creation/division
- Time to die is inherited

Time units are hours.
"""

using Cyton
import Cyton: shouldDie, shouldDivide, inherit, step

using DataFrames
using Gadfly: plot, layer, cm, Gadfly, Theme, Guide, Geom, Col

# Gadfly defaults
Gadfly.set_default_plot_size(20cm, 20cm)
Gadfly.push_theme(Theme(background_color="white"))

# Parameters from the Cyton2 paper
λ_firstDivision = LogNormalParms(log(39.89), 0.28)
λ_subsequentDivision = FixedDistributionParms(9.21)
λ_divisionDestiny = LogNormalParms(log(71.86), 0.11)
λ_lifetime = LogNormalParms(log(116.8), 0.85)

@enum Fate Undivided Dividing Destiny


"""
Time to die is drawn from a distribution when the cell is created.
Daughter cells will inherit this. It is nulled if the cell becomes
a memory cell.
"""
struct DeathTimer <: FateTimer
  timeToDeath::Float64
end
function DeathTimer(r::DistributionParmSet)
  DeathTimer(draw(r))
end
inherit(timer::DeathTimer, time::Float64) = timer
step(timer::DeathTimer, time::Float64, Δt::Float64) = nothing

"Time to divide drawn from distribution"
struct DivisionTimer <: FateTimer
  timeToDivision::Float64
  timeToDestiny::Float64
end
"Constructor for fresh cells"
DivisionTimer(division::DistributionParmSet, destiny::DistributionParmSet) = DivisionTimer(draw(division), draw(destiny))
"Constructor for daughter cells"
DivisionTimer(r::DistributionParmSet, start::Float64, destiny::Float64) = DivisionTimer(draw(r) + start, destiny)
step(timer::DivisionTimer, time::Float64, Δt::Float64) = nothing
inherit(timer::DivisionTimer, time::Float64) = DivisionTimer(λ_subsequentDivision, time, timer.timeToDestiny)

"Create a new cell"
function cellFactory(birth::Float64=0.0)
  cell = Cell(birth)
  addTimer(cell, DeathTimer(λ_lifetime))
  addTimer(cell, DivisionTimer(λ_firstDivision, λ_divisionDestiny))
  return cell
end

"Indicate whether the cell should die and be removed."
shouldDie(death::DeathTimer, time::Float64) = time > death.timeToDeath

"Indicate the cell will divide. Must be earlier than destiny and after the next division time"
shouldDivide(division::DivisionTimer, time::Float64) = error("is this working?")#time < division.timeToDestiny && time > division.timeToDivision

function run(model::CellPopulation, runDuration::Float64)
  print("Time to run:")
  @time begin
    counts = DataFrame(time=Float64[], 
    total=[], 
    gen0 = [],
    gen1 = [],
    gen2 = [],
    gen3 = [],
    gen4 = [],
    gen5 = [],
    gen6 = [],
    gen7 = [],
    gen8 = [],
    genOther = []
    )
    Δt = modelTimeStep(model)
    for tm in 1:Δt:runDuration
      step(model)

      local genCnts = zeros(10)
      cells = model.cells
      for cell in cells
        gen = cell.generation
        if gen <= 8
          genCnts[gen+1] += 1
        else
          genCnts[10] += 1
        end
      end
      push!(counts, (tm, length(cells), genCnts...))
    end
  end

  vars = [:total :gen0 :gen1 :gen2 :gen3 :gen4 :gen5 :gen6 :gen7 :gen8 :genOther]
  h = plot(counts, x=:time, y=Col.value(vars...), color=Col.index(vars...))
  display(h)

  println("Done at model time=$(modelTime(model))")
end

println(rpad(lpad(" start ", 30, "-"), 55, "-"))
model = createPopulation(1000, cellFactory)
run(model, 100.0)
