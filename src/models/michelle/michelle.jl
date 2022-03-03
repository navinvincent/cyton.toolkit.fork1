"""
A cell model based on a mix of Michelle's experiments and the Cyton2 paper.
- Time to first division and destiny are based on the Cyton2 data
- Death is based on a model of protein levels over time. An ensemble value
  is calculated and if the ensemble is less than a threshold the cell dies.
  The ensemble is a weighted sum of the individual protein levels. 

"""

using Cyton
import Cyton: shouldDie, shouldDivide, inherit, step

using DataFrames, Colors, Cairo
using Gadfly: plot, layer, cm, Gadfly, Theme, Guide, Geom, Col, mm, style, Scale

using Base.Threads, Serialization
import Base.show

abstract type Parameters end

# Gadfly defaults
Gadfly.set_default_plot_size(20cm, 20cm)
Gadfly.push_theme(Theme(background_color="white"))

# Parameters from the Cyton2 paper
λ_firstDivision = LogNormalParms(log(39.89), 0.28)
λ_subsequentDivision = FixedDistributionParms(9.21)
λ_divisionDestiny = LogNormalParms(log(71.86), 0.11)
# λ_lifetime = LogNormalParms(log(116.8), 0.85)

# Parameters for the distributions protein level time courses
include("timecourses.jl")
Γ_timeCourseParms(gstd::Float64) = GammaTimeCourseParms(LogNormal(0, gstd), 114.2, 4.6, 24.1)

proteins = ["BCL2", "BCLxL", "MCL1", "BIM"]
initialWeights = Dict{String, Float64}([
  "BCL2"  =>  5,
  "BCLxL" =>  5,
  "MCL1"  =>  5,
  "BIM"   => -5,
])
initialThreshold = 6.0
gracePeriod = 72.0 # hours

#----------------- Death machinery ------------------
struct ThresholdDeath <: FateTimer
  "The threshold below which the cell will die"
  threshold::Float64
  "The protein level weights"
  weights::Dict{String, Float64}
  "References to the current protein level"
  proteinLevels::Dict{String, TimeCourseParameters}
  "Period before death kicks in"
  gracePeriod::Float64
end
function ThresholdDeath(threshold::Float64, weights::Dict{String, Float64}, gracePeriod::Float64, levelGenerator::TimeCourseParameters)
  pl = Dict{String, TimeCourseParameters}()
  for p in proteins
    pl[p] = levelGenerator
  end
  ThresholdDeath(threshold, weights, pl, gracePeriod)
end

ensemble(::FateTimer, ::Float64) = 0.0
function ensemble(death::ThresholdDeath, time::Float64) 
  e = 0.0
  for (p, w) in death.weights
    level = death.proteinLevels[p]
    e += w * level(time)
  end
  return e
end
function ensemble(cell::Cell, time::Float64)
  e = 0.0
  for timer in cell.timers
    e += ensemble(timer, time)
  end
  return e
end

function proteinLevel(cell::Cell, protein::String, time::Float64)
  l = 0.0
  for timer in cell.timers
    l += proteinLevel(timer, protein, time)
  end
  return l
end
proteinLevel(::FateTimer, ::String, ::Float64) = 0.0
proteinLevel(death::ThresholdDeath, protein::String, time::Float64) = death.proteinLevels[protein](time)

function shouldDie(death::ThresholdDeath, time::Float64)
  if time < death.gracePeriod
    return false
  end

  return ensemble(death, time) < death.threshold
end

step(timer::ThresholdDeath, time::Float64, Δt::Float64) = nothing

"Daughters inherit the protein levels??"
inherit(deathTimer::ThresholdDeath, time::Float64) = deathTimer
#-----------------------------------------------------

#--------------- Division machinery ------------------
"Time to divide drawn from distribution"
struct DivisionTimer <: FateTimer
  timeToDivision::Float64
  timeToDestiny::Float64
end

"Constructor for fresh cells"
DivisionTimer(division::DistributionParmSet, destiny::DistributionParmSet) = DivisionTimer(draw(division), draw(destiny))

step(timer::DivisionTimer, time::Float64, Δt::Float64) = nothing

"Daughter cells get a new division timer and inherit the destiny timer"
inherit(timer::DivisionTimer, time::Float64) = DivisionTimer(λ_subsequentDivision, time, timer.timeToDestiny)
DivisionTimer(r::DistributionParmSet, start::Float64, destiny::Float64) = DivisionTimer(draw(r) + start, destiny)

"Indicate the cell will divide. Must be earlier than destiny and after the next division time"
shouldDivide(division::DivisionTimer, time::Float64) = time < division.timeToDestiny && time > division.timeToDivision
#----------------------------------------------------

#------------------ Cell factory --------------------
function cellFactory(parameters::Parameters, birth::Float64=0.0)
  michelleCell = Cell(birth)

  death = ThresholdDeath(parameters.threshold, initialWeights, gracePeriod, Γ_timeCourseParms(parameters.gstd))
  addTimer(michelleCell, death)
  division = DivisionTimer(λ_firstDivision, λ_divisionDestiny)
  addTimer(michelleCell, division)

  return michelleCell
end
#----------------------------------------------------

struct Result
  counts::DataFrame
  proteinLevels::DataFrame
end

function run(model::CellPopulation, runDuration::Float64)
  Δt = modelTimeStep(model)
  time = 0:Δt:runDuration
  count = zeros(Int, length(time))

  proteinSampleTimes = Set([72.0, 96.0, 120.0, 144.0, 168.0])
  proteinLevels = DataFrame(time=Float64[], name=String[], level=Float64[])
  for (i, tm) in enumerate(time)
    step(model)
    count[i] = cellCount(model)

    if tm % 10 == 0
      println("$(tm): $(count[i]) cells")
    end

    if tm in proteinSampleTimes
      
      for protein in proteins
        for cell in model.cells
          level = proteinLevel(cell, protein, tm)
          push!(proteinLevels, (tm, protein, level))
        end
      end

      for cell in model.cells
        level = ensemble(cell, tm)
        push!(proteinLevels, (tm, "ensemble", level))
      end
    end
  end

  counts = DataFrame(time=time, count=count)
  result = Result(counts, proteinLevels)

  return result
end

struct ConcreteParameters <: Parameters
  threshold::Float64
  gstd::Float64
end
show(io::IO, ::Parameters) = print(io, "no desciption")
show(io::IO, parms::ConcreteParameters) = print(io, "threshold=$(parms.threshold) gstd=$(parms.gstd)")

lk = ReentrantLock()
results = Dict{Parameters, Result}()

thresholds = 5.0:8.0
gstds      = 0.1:0.1:0.5
parameters = ConcreteParameters[]
for threshold in thresholds
  for gstd in gstds
    p = ConcreteParameters(threshold, gstd)
    push!(parameters, p)
  end
end

results = Dict{ConcreteParameters, Result}()
@threads for parameter in parameters
  model = createPopulation(20000, (birth) -> cellFactory(parameter, birth))
  result = run(model, 200.0);
  lock(lk) do 
    results[parameter] = result
  end
  @info("$(parameter) done")
end

@info("Runs finished!")

open("results.dat", "w") do io
  serialize(io, results)
end

@info("Data saved!")

include("plotting.jl")

@info("Data plotted!")

@info("You are awesome!")
