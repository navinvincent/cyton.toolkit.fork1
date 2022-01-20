
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
  shouldDivide,
  Cell,
  draw


"A mode of around 24 hours"
λ_firstDivision = LogNormalParms(3.2, 0.4)
"A mode of around 6 hours"
λ_subsequentDivision = LogNormalParms(2.0, 0.4)
"A mode of around 48 hours"
λ_lifetime = LogNormalParms(4.0, 0.4)
"A mode of around 5 days"
λ_divisionDestiny = LogNormalParms(5.0, 0.5)
"A constant rate of conversion to memory cells"
memoryCellRate = 0.001 # Conversions per hour

@enum CellType PreFirstDivision Dividing Destiny Memory

"""
This differentiation accumulator does everything. It holds the conversion to memory cells
and the times to differentiate from pre division -> dividing -> destiny. It therefore needs
to track the cell type. An alternative mechanism could be to use subtypes of this class to
map the cell types.
"""
struct ConstantDifferentiator <: DifferentationAccumulator
  "Time of next differentiation event"
  nextEvent::Float64
  "Constant rate of conversion to memory cells"
  memoryCellRate::Float64
  "Cell will convert to memory on next time step"
  convertToMemory::Bool
  "The cell type"
  cellType::CellType

  function ConstantDifferentiator(r::DistributionParmSet, cellType::CellType)
    new(draw(r), memoryCellRate, false, cellType)
  end
end

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
struct SimpleDeath <: DeathAccumulator
  timeToDie::Float64
  function SimpleDeath(r::DistributionParmSet)
    new(draw(r))
  end
end

"Time to divide drawn from distribution"
struct SimpleDivision <: DivisionAccumulator
  timeToDivision::Float64
  function SimpleDivision(r::DistributionParmSet, startTime::Float64)
    new(draw(r) + startTime)
  end
end

"Create a new cell"
function cellFactory(birth::Float64=0.0)
  cell = Cell(birth)
  cell.differentiationAccumulator = ConstantDifferentiator(λ_firstDivision, PreFirstDivision)
  cell.deathAccumulator = SimpleDeath(λ_lifetime)
  cell.divisionAccumulator = nothing
  return cell
end

"Indicate whether the cell should die and be removed."
shouldDie(cell::Cell, time::Float64) = age(cell, time) > cell.deathAccumulator.timeToDie 

"Indicate the cell will divide."
shouldDivide(divisionA::SimpleDivision, time::Float64) = time > divisionA.timeToDivision

"Indicate the cell will differentiate"
function shouldDifferentiate(differentiateA::ConstantDifferentiator, time::Float64)
  if differentiateA.convertToMemory
    return true
  end

  return differentiateA.nextEvent > time
end

"""
Enact the differentiation state machine
Any cell can potentially becomes a memory cell.
Other transitions are:
pre dividing -> dividing
dividing -> destiny
"""
function differentiate(cell::Cell, time::Float64)
  if cell.differentiateA.convertToMemory
    cell.deathAccumulator = nothing
    cell.divisionAccumulator = nothing
    cell.differentiationAccumulator.cellType = Memory
    return
  end

  cellType = cell.differentiationAccumulator.cellType

  if cellType == PreFirstDivision
    cell.divisionAccumulator = SimpleDivision(λ_subsequentDivision, time)
    cell.differentiationAccumulator.cellType = Dividing
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
  new_cell.divisionAccucumulator = SimpleDivision(λ_subsequentDivision, time)

  return new_cell
end
