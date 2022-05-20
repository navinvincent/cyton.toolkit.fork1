
export FateTimer, addTimer, Cell, Stimulus, stimulate, age, CellEvent, Division, Death, addObserver, CellType, cellType, GenericCell

#---------------------- Cell events -----------------------
"""
Top level type for modelling cellular events. 
`Death` and `Division` are handled internally but users should add additional
events to model additional behaviour.
"""
abstract type CellEvent end
"A timer (or something else) has triggered cell division"
struct Division <: CellEvent end
"A timer (or something else) has triggered a cell death event"
struct Death <: CellEvent end
#----------------------------------------------------------

#--------------------- Cell data types --------------------
"Top level Cell type"
abstract type CellType end
"A basic cell type, useful where there is no cell type differentiation in the model"
struct GenericCell <: CellType end
"This type associates another type with each cell, e.g. phenotype, genotype, etc"
abstract type AbstractCell{T <: CellType} end

"""
The main user facing Cell data type. It encapsulates another 
type, `T`, which can model concepts such as phenotype or genotype
"""
mutable struct Cell{T} <: AbstractCell{T}
  "The time the cell is created."
  birth::Float64
  "The generation number of this cell."
  generation::Int64
  "Timers determing the fate of this cell"
  timers::Vector{FateTimer}
  "Observers watching for cell events"
  observers::Dict{CellEvent, Vector{Function}}
  "An arbitrary type "
  cellType::T
end

# Convenience constructors
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
#----------------------------------------------------------

#------------------------ Stimuli -------------------------
"Top level type for modelling cell stimuli"
abstract type Stimulus end
#----------------------------------------------------------

#-------- Types for wrapping the Agent framework ----------
mutable struct CellAgent <: AbstractAgent
  "These are required by the ABM framework"
  id::Int
  pos::NTuple{2, Int}
  cell::AbstractCell
end
CellAgent(cell::AbstractCell) = CellAgent(0, (0, 0), cell)
CellAgent(id::Int, cell::AbstractCell) = CellAgent(id, (0, 0), cell)
#----------------------------------------------------------

function stimulate(::Cell, ::Stimulus, time::Time, Δt::Duration)::Nothing end
addTimer(cell::Cell, timer::FateTimer) = push!(cell.timers, timer)

age(cell::Cell, time::Time) = time - cell.birth

"Internal, probably not necassary"
function die(::AbstractCell)::Nothing end

"""
Create a daughter from the mother and reset the mother's state.
This is called internally.
"""
function divide(cell::AbstractCell, time::Time) 
  # Reset the state of the mother
  cell.generation += 1
  cell.birth = time

  new_cell = Cell(time, cell.generation)

  # new and existing cell "inherit" the mother's timers 
  # see the inherit method for more detail
  for i in 1:length(cell.timers)
    timer = cell.timers[i]
    cell.timers[i] = inherit(timer, time)
    addTimer(new_cell, inherit(timer, time))
  end

  return new_cell
end

#----------------- observers --------------------
"""
The adds a callback function,
  `observer(event::CellEvent, cell::Cell, time::Time)`, 
to a cell. This function is called when `Cell` triggers a `CellEvent`.
"""
function addObserver(event::CellEvent, cell::Cell, observer::Function) 
  observers = cell.observers
  if !haskey(observers, event)
    observers[event] = Function[]
  end
  push!(observers[event], observer)
end

"This is an internal function."
function notifyObservers(event::CellEvent, cell::Cell, time::Time)
  for observer in get(cell.observers, event, Function[])
    observer(event, cell, time)
  end
end
#-----------------------------------------------
