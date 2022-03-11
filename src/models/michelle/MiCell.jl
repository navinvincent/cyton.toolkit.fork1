
using Cyton

import Cyton: shouldDie, shouldDivide, inherit, step

include("timecourses.jl")


abstract type Parameters end

# Parameters from the Cyton2 paper
λ_firstDivision = LogNormalParms(log(39.89), 0.28)
λ_subsequentDivision = FixedDistributionParms(9.21)
λ_divisionDestiny = LogNormalParms(log(71.86), 0.11)
# λ_lifetime = LogNormalParms(log(116.8), 0.85)

# Parameters for the distributions protein level time courses
include("timecourses.jl")
# TimeCourseParms(gstd::Float64) = GammaTimeCourseParms(LogNormal(0, gstd), 114.2, 4.6, 24.1)
TimeCourseParms(gstd::Float64) = PiecewiseLinear(LogNormal(0, gstd), 235.0)

proteins = ["BCL2", "BCLxL", "MCL1", "BIM"]
initialWeights = Dict{String, Float64}([
  "BCL2"  =>  5,
  "BCLxL" =>  5,
  "MCL1"  =>  5,
  "BIM"   => -5,
])
initialThreshold = 6.0
gracePeriod = 72.0

#----------------- Death machinery ------------------
mutable struct ThresholdDeath <: FateTimer
  "The threshold below which the cell will die"
  threshold::Float64
  "The protein level weights"
  weights::Dict{String, Float64}
  "References to the current protein level"
  proteinLevels::Dict{String, TimeCourseParameters}
  "A toggle to indicate the cell can now die"
  canDie::Bool
  "A mechanismn to keep using the timer even after the cell nominally died"
  isDead::Bool
end
function ThresholdDeath(threshold::Float64, weights::Dict{String, Float64}, levelGenerator::TimeCourseParameters)
  pl = Dict{String, TimeCourseParameters}()
  for p in proteins
    pl[p] = levelGenerator
  end
  ThresholdDeath(threshold, weights, pl, false, false)
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
  # global gracePeriod

  # if time < gracePeriod
  #   return false
  # end

  if death.isDead
    return true
  end

  e = ensemble(death, time)
  t = death.threshold
  if e > t
    death.canDie = true
  end

  if death.canDie && e < t
    death.isDead = true
    return true
  else
    return false
  end
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
