
export FateTimer, addTimer, Cell, Stimulus, stimulate, age, die, CellEvent, Division, Death, addObserver, CellType, cellType

abstract type CellType end
struct GenericCell <: CellType end
abstract type AbstractCell{T <: CellType} end
abstract type Stimulus end

abstract type CellEvent end
struct Division <: CellEvent end
struct Death <: CellEvent end

mutable struct CellAgent <: AbstractAgent
  "These are required by the ABM framework"
  id::Int
  pos::NTuple{2, Int}
  cell::AbstractCell
end
CellAgent(cell::AbstractCell) = CellAgent(0, (0, 0), cell)
CellAgent(id::Int, cell::AbstractCell) = CellAgent(id, (0, 0), cell)

mutable struct Cell{T} <: AbstractCell{T}
  birth::Float64
  generation::Int64
  timers::Vector{FateTimer}
  observers::Dict{CellEvent, Vector{Function}}
  cellType::T
end

function Cell(birth::Float64)
  return Cell(birth, 0, FateTimer[], Dict{CellEvent, Vector{Function}}(), GenericCell())
end

function Cell(birth::Float64, cellType::T) where T <: CellType
  return Cell(birth, 0, FateTimer[], Dict{CellEvent, Vector{Function}}(), cellType)
end

function Cell(birth::Float64, divisionCount::Int64)
  return Cell(birth, divisionCount, FateTimer[], Dict{CellEvent, Vector{Function}}(), GenericCell())
end

cellType(::AbstractCell{T}) where T <: CellType = T

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

#----------------- observers --------------------
function addObserver(event::CellEvent, cell::Cell, observer::Function) 
  observers = cell.observers
  if !haskey(observers, event)
    observers[event] = Function[]
  end
  push!(observers[event], observer)
end

function notifyObservers(event::CellEvent, cell::Cell, time::Float64)
  for observer in get(cell.observers, event, Function[])
    observer(event, cell, time)
  end
end
#-----------------------------------------------
