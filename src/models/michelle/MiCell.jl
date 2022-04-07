
using Cyton

import Cyton: shouldDie, shouldDivide, inherit, step

include("utils.jl")

struct WildType <: CellType end
struct WildTypeDrugged <: CellType end
struct KO <: CellType end

# Parameters from the Cyton2 paper
λ_firstDivision = LogNormalParms(log(39.89), 0.28)
λ_subsequentDivision = FixedDistributionParms(9.21)
λ_divisionDestiny = LogNormalParms(log(71.86), 0.11)
# λ_lifetime = LogNormalParms(log(116.8), 0.85)

# Parameters for the distributions protein level time courses
TimeCourseParms(gstd::Float64) = GammaTimeCourseParms(LogNormal(0, gstd), 114.2, 4.6, 24.1)
# TimeCourseParms(gstd::Float64) = PiecewiseLinear(LogNormal(0, gstd), 235.0)

proteins = ["BCL2", "BCLxL", "MCL1", "BIM"]
initialWeights = Dict{String, Float64}([
  "BCL2"  =>  5,
  "BCLxL" =>  1,
  "MCL1"  =>  5,
  "BIM"   => -5,
])
initialThreshold = 6.0
gracePeriod = 0.0

#----------------- Death machinery ------------------
mutable struct ThresholdDeath <: FateTimer
  "The threshold below which the cell will die"
  threshold::Float64
  "The protein level weights"
  weights::Dict{String, Float64}
  "References to the current protein level"
  proteinLevels::Dict{String, TimeCourseParameters}
  "The cell cannot die before this time"
  gracePeriod::Float64
  "Toggle to allow death after first division"
  canDie::Bool
end
function ThresholdDeath(threshold::Float64, weights::Dict{String, Float64}, levelGenerator::Function)
  pl = Dict{String, TimeCourseParameters}()
  for p in proteins
    pl[p] = levelGenerator()
  end
  ThresholdDeath(threshold, weights, pl, gracePeriod, false)
end

ensemble(::FateTimer, ::Float64) = 0.0
function ensemble(death::ThresholdDeath, time::Float64) 
  e = 0.0
  for (p, w) in death.weights
    level = death.proteinLevels[p]
    e += w * level(time)
  end
  return maximum([0.0, e])
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

  e = ensemble(death, time)
  t = death.threshold

  return e < t
end

step(timer::ThresholdDeath, time::Float64, Δt::Float64) = nothing

"Daughters inherit the protein levels"
function inherit(deathTimer::ThresholdDeath, time::Float64) 
  deathTimer.canDie = true
  deathTimer
end

isDead(timer::FateTimer) = false
isDead(timer::ThresholdDeath) = timer.isDead
function isDead(cell::Cell)
  for timer in cell.timers
    if isDead(timer)
      return true
    end
  end
  return false
end
#-----------------------------------------------------

#---------------- Division machinery ------------------
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
#------------------------------------------------------


#---------------- BH3 mimetic treatment ---------------
struct Bh3Stimulus <: Stimulus
  applicationTime::Float64
  inhibitionFactor::Float64
  protein::String
end

function stimulate(timer::FateTimer, stimulus::Stimulus, time::Float64) end
function stimulate(death::ThresholdDeath, bh3::Bh3Stimulus, time::Float64)
  if time == bh3.applicationTime
    death.weights[bh3.protein] *= bh3.inhibitionFactor
  end
end

function stimulate(cell::Cell{WildTypeDrugged}, bh3::Bh3Stimulus, time::Float64)
  for timer in cell.timers
    stimulate(timer, bh3, time)
  end
end
#------------------------------------------------------

