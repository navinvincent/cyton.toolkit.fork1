module Cyton

using Agents

export modelTime, modelTimeStep, step, createPopulation, CellPopulation, cellCount

include("probability/DistributionParms.jl")
include("cells/CellModules.jl")
include("cells/Cells.jl")
include("utils/void_space.jl")

noopCallback(::Cell, ::Float64) = nothing
mutable struct CellPopulation
  model::AgentBasedModel
  deathCallback::Function
  divisionCallback::Function
end
CellPopulation(model::AgentBasedModel) = CellPopulation(model, noopCallback, noopCallback)

function Base.getproperty(population::CellPopulation, v::Symbol)
  if v == :cells
    return Iterators.map(a -> a.cell, values(population.model.agents))
  end
  return getfield(population, v) # Just fall through for other fields
end

"Get the current model time"
modelTime(model::CellPopulation) = model.model.properties[:step_cnt] * modelTimeStep(model)
"Get the model time step"
modelTimeStep(model::CellPopulation) = model.model.properties[:Δt]
"Get the number of cells in the population"
cellCount(model::CellPopulation) = length(model.model.agents)

"""
Create a population of cells:
 nCells: size of starting populations
 cellFactory: A function that retuns a cell
 deathCallback: A function that is called just before a cell dies
 divisionCallback: A function that is call just before a cell divides
"""
function createPopulation(nCells::Int, 
  cellFactory::Function; 
  deathCallback::Function=noopCallback,
  divisionCallback::Function=noopCallback)

  space = VoidSpace()
  scheduler = Schedulers.fastest
  properties = Dict(:step_cnt => 0, :Δt => 0.1)
  model = AgentBasedModel(CellAgent, space; properties, scheduler)

  for id in 1:nCells
    cell = cellFactory(0.0)
    agent = CellAgent(id, cell)
    add_agent_pos!(agent, model)
  end

  return CellPopulation(model, deathCallback, divisionCallback)
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

function doStep(agent::CellAgent, time::Float64, Δt::Float64, model::CellPopulation, stimuli::Vector{T}) where T<:Stimulus
  cell = agent.cell

  willDie = false
  willDivide = false
  for timer in cell.timers
    step(timer, time, Δt)
    willDie = willDie || shouldDie(timer, time) 
    willDivide = willDivide || shouldDivide(timer, time) 
  end

  if willDie
    model.deathCallback(cell, time)
    die(cell)
    kill_agent!(agent, model.model)
  end

  for stimulus in stimuli
    stimulate(cell, stimulus, time)
  end

  if willDivide
    model.divisionCallback(cell, time)
    new_cell = divide(cell, time)
    if new_cell ≠ nothing
      new_agent = CellAgent(model.model.maxid[]+1, new_cell)
      add_agent_pos!(new_agent, model.model)
    end
  end
end

end
