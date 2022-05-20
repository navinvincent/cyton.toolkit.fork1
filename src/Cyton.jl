module Cyton

using Agents
import Base.length

export modelTime, modelTimeStep, step, createPopulation, cellCount, CellPopulation, cohortCount, Time, Duration

const Time = Float64
const Duration = Float64

include("probability/DistributionParms.jl")
include("cells/CellModules.jl")
include("cells/Cells.jl")
include("utils/void_space.jl")

#------------------------ Cell population -----------------------
"""
The population of cells with convenince constructors. Ordinarily, this will be constructed
by the framework.
"""
mutable struct CellPopulation
  model::AgentBasedModel
  startingCnt::Int
  eventCallbacks::Vector{Function} 
end
CellPopulation(model::AgentBasedModel) = CellPopulation(model, length(model.agents), Function[])
CellPopulation(model::AgentBasedModel, eventCallbacks::Vector{Function}) = 
CellPopulation(model, length(model.agents), eventCallbacks)


function Base.getproperty(population::CellPopulation, v::Symbol)
  if v == :cells
    return Iterators.map(a -> a.cell, values(population.model.agents))
  end
  return getfield(population, v) # Just fall through for other fields
end

Base.length(population::CellPopulation) = length(population.model.agents)
#----------------------------------------------------------------

"Get the current model time"
modelTime(model::CellPopulation) = model.model.properties[:step_cnt] * modelTimeStep(model)
"Get the model time step"
modelTimeStep(model::CellPopulation) = model.model.properties[:Δt]
"Get the number of cells in the population"
cellCount(model::CellPopulation) = length(model.model.agents)
"The current cohort count, cell count normalised by generation number"
function cohortCount(model::CellPopulation)
  cohort = 0.0
  for cell in model.cells
    cohort += 2.0^-cell.generation
  end
  return cohort/model.startingCnt
end

"""
Create a population of cells:
 nCells: size of starting populations
 cellFactory: A function that returns constructs a new cell
 eventCallbacks: Function that are called when events occurs
 
 `Division` and `Death` are two predefined events
"""
function createPopulation(nCells::Int, 
  cellFactory::Function; 
  eventCallbacks::Vector{Function}=Function[])

  space = VoidSpace()
  scheduler = Schedulers.fastest
  properties = Dict(:step_cnt => 0, :Δt => 0.1)
  model = AgentBasedModel(CellAgent, space; properties, scheduler)

  for id in 1:nCells
    cell = cellFactory(0.0)
    agent = CellAgent(id, cell)
    add_agent_pos!(agent, model)
  end

  return CellPopulation(model, eventCallbacks)
end

"Step the population forward in time by one time step, with stimulus"
step(model::CellPopulation, stimulus::T) where T<:Stimulus = step(model, [stimulus])
"Step the population forward in time by one time step, with optional stimuli"
step(model::CellPopulation, stimuli::Vector{T}=Vector{Stimulus}()) where T<:Stimulus = step!(model.model, (a, _) -> step(a, model, stimuli), stepModel)

"Step a cell forward in time by one time step"
function step(agent::CellAgent, model::CellPopulation, stimuli::Vector{T}) where T<:Stimulus
  Δt   = modelTimeStep(model)
  time = modelTime(model)
  doStep(agent, time, Δt, model, stimuli)
end

stepModel(model::AgentBasedModel) = model.properties[:step_cnt] += 1

function doStep(agent::CellAgent, time::Time, Δt::Duration, model::CellPopulation, stimuli::Vector{T}) where T<:Stimulus
  cell = agent.cell

  for stimulus in stimuli
    stimulate(cell, stimulus, time, Δt)
  end

  events = map(cell.timers) do timer
    step(timer, time, Δt)
  end
  events = filter(x -> x ≠ nothing, events)

  if any(typeof.(events) .== Death)
    die(cell)
    kill_agent!(agent, model.model)
  end
  
  if any(typeof.(events) .== Division)
    new_cell = divide(cell, time)
    if new_cell ≠ nothing
      new_agent = CellAgent(model.model.maxid[]+1, new_cell)
      add_agent_pos!(new_agent, model.model)
      for e in events
        notifyObservers(e, new_cell, time)
      end
    end
  end

  for e in events
    # Cell observers
    notifyObservers(e, cell, time)
    # Population observers
    for cb in model.eventCallbacks
      cb(e, time)
    end
  end

end

end
