
abstract type AbstractCell end
abstract type Stimulus end

export FateTimer, addTimer, Cell, Stimulus, stimulate, age, die

mutable struct CellAgent <: AbstractAgent
  "These are required by the ABM framework"
  id::Int
  pos::NTuple{2, Int}
  cell::AbstractCell
end
CellAgent(cell::AbstractCell) = CellAgent(0, (0, 0), cell)
CellAgent(id::Int, cell::AbstractCell) = CellAgent(id, (0, 0), cell)

mutable struct Cell <: AbstractCell
  birth::Float64
  generation::Int64
  timers::Vector{FateTimer}
end
  
function Cell(birth::Float64)
  return Cell(birth, 0, [])
end

function Cell(birth::Float64, divisionCount::Int64)
  return Cell(birth, divisionCount, [])
end

stimulate(_::Cell, _::Stimulus, time::Float64) = nothing
addTimer(cell::Cell, timer::FateTimer) = push!(cell.timers, timer)

age(cell::Cell, time::Float64) = time - cell.birth

die(cell::AbstractCell) = nothing

"Create a daughter from the mother and reset the mother's state"
function divide(cell::AbstractCell, time::Float64) 
  cell.generation += 1
  new_cell = Cell(time, cell.generation)

  for i in 1:length(cell.timers)
    timer = cell.timers[i]
    cell.timers[i] = inherit(timer, time)
    addTimer(new_cell, inherit(timer, time))
  end

  return new_cell
end

