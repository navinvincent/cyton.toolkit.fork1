
using Cyton

import Cyton: shouldDie, shouldDivide, inherit, step

include("utils.jl")

struct WildType <: CellType end
struct WildTypeDrugged <: CellType end
struct KO <: CellType end

# Parameters from the Cyton2 paper
# λ_firstDivision = LogNormalParms(log(39.89), 0.28)
# λ_subsequentDivision = FixedDistributionParms(9.21)
# λ_divisionDestiny = LogNormalParms(log(71.86), 0.11)
# λ_lifetime = LogNormalParms(log(116.8), 0.85)
# Parameters from fitting MR-70 with cyton solver
λ_firstDivision = LogNormalParms(log(46.2), 0.37)
λ_subsequentDivision = FixedDistributionParms(11.1)
λ_divisionDestiny = LogNormalParms(log(74.4), 0.08)

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
gracePeriod = 72

#----------------- Death machinery ------------------
mutable struct ThresholdDeath <: FateTimer
  "The threshold below which the cell will die"
  threshold::Float64
  "The protein level weights"
  weights::Dict{String, Float64}
  "References to the current protein level"
  proteinLevels::Dict{String, TimeCourseParameters}
  "The cell cannot die before this time"
  gracePeriod::Duration
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
function ensemble(death::ThresholdDeath, time::Time) 
  e = 0.0
  for (p, w) in death.weights
    level = death.proteinLevels[p]
    e += w * level(time)
  end
  return maximum([0.0, e])
end
function ensemble(cell::Cell, time::Time)
  e = 0.0
  for timer in cell.timers
    e += ensemble(timer, time)
  end
  return e
end

function proteinLevel(cell::Cell, protein::String, time::Time)
  l = 0.0
  for timer in cell.timers
    l += proteinLevel(timer, protein, time)
  end
  return l
end
proteinLevel(::FateTimer, ::String, ::Time) = 0.0
proteinLevel(death::ThresholdDeath, protein::String, time::Time) = death.proteinLevels[protein](time)

function shouldDie(death::ThresholdDeath, time::Time)
  if time < death.gracePeriod
    return false
  end

  e = ensemble(death, time)
  t = death.threshold

  return e < t
end

function step(timer::ThresholdDeath, time::Time, Δt::Duration) 
  if shouldDie(timer, time)
    return Death()
  else
    return nothing
  end
end

"Daughters inherit the protein levels"
function inherit(deathTimer::ThresholdDeath, time::Time) 
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
  timeToDivision::Time
  timeToDestiny::Time
end

"Constructor for fresh cells"
DivisionTimer(division::DistributionParmSet, destiny::DistributionParmSet) = DivisionTimer(draw(division), draw(destiny))

function step(timer::DivisionTimer, time::Time, ::Duration) 
  if shouldDivide(timer, time)
    return Division()
  else
    return nothing
  end
end

"Daughter cells get a new division timer and inherit the destiny timer"
inherit(timer::DivisionTimer, time::Time) = DivisionTimer(λ_subsequentDivision, time, timer.timeToDestiny)
DivisionTimer(r::DistributionParmSet, start::Time, destiny::Float64) = DivisionTimer(draw(r) + start, destiny)

"Indicate the cell will divide. Must be earlier than destiny and after the next division time"
shouldDivide(division::DivisionTimer, time::Time) = time < division.timeToDestiny && time > division.timeToDivision
#------------------------------------------------------


#---------------- BH3 mimetic treatment ---------------
struct Bh3Stimulus <: Stimulus
  applicationTime::Time
  inhibitionFactor::Float64
  protein::String
end

function stimulate(timer::FateTimer, stimulus::Stimulus, time::Time) end
function stimulate(death::ThresholdDeath, bh3::Bh3Stimulus, time::Time)
  if time == bh3.applicationTime
    death.weights[bh3.protein] *= bh3.inhibitionFactor
  end
end

function stimulate(cell::Cell{WildTypeDrugged}, bh3::Bh3Stimulus, time::Time, Δt::Duration)
  for timer in cell.timers
    stimulate(timer, bh3, time)
  end
end
#------------------------------------------------------

