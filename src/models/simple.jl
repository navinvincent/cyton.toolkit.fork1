
"""
A simple cell has the following behaviour
- It does not divide prior to stimulation
- After stimulation there is a delay to first division
- After first division the cell will divide until it reaches destiny
- Most cells will die after a period of time.
- A proportion of cells will not die - memory cells
  - differentiation to a memory cell occurs at constant (slow) rate
  - memory cells don't divide

In this simple, proof of concept model:
- all lifetimes are drawn from LogNormalParms
- Time to next division are initiated de novo on creation/division
- Time to die is inherited
- 

Time units are hours.
"""

import Cyton: 
  DistributionParmSet,
  LogNormalParms,
  DeathAccumulator,
  CellModule,
  DifferentationAccumulator,
  DivisionAccumulator,
  DeathAccumulator,
  SimpleDecay,
  shouldDie,
  shouldDifferentiate,
  differentiate,
  shouldDivide,
  divide,
  Cell,
  draw,
  step

import Base: copy

# Parameters from the Cyton2 paper, transformed according to
# https://discourse.julialang.org/t/lognormal-distribution-how-to-set-mu-and-sigma/7101

λ_firstDivision = LogNormalParms(log(38.4)-0.13^2/2, 0.13)
λ_subsequentDivision = LogNormalParms(log(10.9)-0.23^2/2, 0.23)
λ_divisionDestiny = LogNormalParms(log(53.79)-0.21^2/2, 0.21)
λ_lifetime = LogNormalParms(log(86.66)-0.18^2/2, 0.18)

"A constant rate of conversion to memory cells"
memoryCellRate = 0.00 # Conversions per hour

@enum CellType Undivided Dividing Destiny Memory

"""
This differentiation accumulator does everything. It holds the conversion to memory cells
and the times to differentiate from pre division -> dividing -> destiny. It therefore needs
to track the cell type. An alternative mechanism could be to use subtypes of this class to
map the cell types.
"""
mutable struct ConstantDifferentiator <: DifferentationAccumulator
  "Time of next differentiation event"
  nextEvent::Float64
  "Constant rate of conversion to memory cells"
  memoryCellRate::Float64
  "Cell will convert to memory on next time step"
  convertToMemory::Bool
  "The cell type"
  cellType::CellType
end
function ConstantDifferentiator(r::DistributionParmSet, cellType::CellType)
  ConstantDifferentiator(draw(r), memoryCellRate, false, cellType)
end
copy(x::ConstantDifferentiator) = ConstantDifferentiator(x.nextEvent, x.memoryCellRate, x.convertToMemory, x.cellType)

"The step function for differentiation"
function step(differentiateA::ConstantDifferentiator, time::Float64, Δt::Float64)
  if rand() < differentiateA.memoryCellRate * Δt
    differentiateA.convertToMemory = true
  end
end

"""
Time to die is drawn from a distribution when the cell is created.
Daughter cells will inherit this. It is nulled if the cell becomes a
memory cell.
"""
struct TimedDeath <: DeathAccumulator
  timeToDie::Float64
end
function TimedDeath(r::DistributionParmSet)
  TimedDeath(draw(r))
end
copy(x::TimedDeath) = x
step(_::TimedDeath, time::Float64, Δt::Float64) = nothing

"Time to divide drawn from distribution"
struct TimedDivision <: DivisionAccumulator
  timeToDivision::Float64
  function TimedDivision(r::DistributionParmSet, startTime::Float64)
    new(draw(r) + startTime)
  end
end
step(_::TimedDivision, time::Float64, Δt::Float64) = nothing

"Create a new cell"
function cellFactory(birth::Float64=0.0)
  cell = Cell(birth)
  cell.differentiationAccumulator = ConstantDifferentiator(λ_firstDivision, Undivided)
  cell.deathAccumulator = TimedDeath(λ_lifetime)
  cell.divisionAccumulator = nothing
  return cell
end

"Indicate whether the cell should die and be removed."
shouldDie(death::TimedDeath, time::Float64) = time > death.timeToDie 

"Indicate the cell will divide."
shouldDivide(divisionA::TimedDivision, time::Float64) = time > divisionA.timeToDivision

"Indicate the cell will differentiate"
function shouldDifferentiate(differentiateA::ConstantDifferentiator, time::Float64)
  if differentiateA.convertToMemory
    return true
  end

  return time > differentiateA.nextEvent
end

"""
Enact the differentiation state machine
Any cell can potentially becomes a memory cell.
Other transitions are:
pre dividing -> dividing
dividing -> destiny
"""
function differentiate(cell::Cell, time::Float64)
  cellType = cell.differentiationAccumulator.cellType

  if cellType == Memory
    return
  end

  if cell.differentiationAccumulator.convertToMemory
    cell.deathAccumulator = nothing
    cell.divisionAccumulator = nothing
    cell.differentiationAccumulator.cellType = Memory
    return
  end


  if cellType == Undivided
    cell.divisionAccumulator = TimedDivision(λ_subsequentDivision, time)
    cell.differentiationAccumulator.cellType = Dividing
    cell.differentiationAccumulator.nextEvent = draw(λ_divisionDestiny)
    return
  end

  if cellType == Dividing
    cell.divisionAccumulator = nothing
    cell.differentiationAccumulator.cellType = Destiny
    return
  end
end

"Create a daughter from the mother and reset the mother's state"
function divide(cell::Cell, time::Float64)
  cell.divisonCount += 1
  cell.birth = time
  new_cell = Cell(cell)
  
  "Inherit from mother"
  new_cell.deathAccumulator = copy(cell.deathAccumulator)
  new_cell.differentiationAccumulator = copy(cell.differentiationAccumulator)

  "de novo timers"
  cell.divisionAccumulator = TimedDivision(λ_subsequentDivision, time)
  new_cell.divisionAccumulator = TimedDivision(λ_subsequentDivision, time)

  return new_cell
end
