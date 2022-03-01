module Cyton

using Agents

export modelTime, modelTimeStep, step, createModel, CellPopulation

include("probability/DistributionParms.jl")
include("cells/CellModules.jl")
include("cells/Cells.jl")
include("utils/void_space.jl")

struct CellPopulation
  model::AgentBasedModel
end
function Base.getproperty(population::CellPopulation, v::Symbol)
  if v == :cells
    return Iterators.map(a -> a.cell, values(population.model.agents))
  end
  return getfield(population, v) # Will probability error
end

modelTime(model::CellPopulation) = model.model.properties[:step_cnt] * modelTimeStep(model)
modelTimeStep(model::CellPopulation) = model.model.properties[:Δt]

function createModel(Ncells::Int, cellFactory::Function)
  space = VoidSpace()
  scheduler = Schedulers.fastest
  properties = Dict(:step_cnt => 0, :Δt => 0.1)
  model = AgentBasedModel(CellAgent, space; properties, scheduler)

  for id in 1:Ncells
    cell = cellFactory()
    agent = CellAgent(id, cell)
    add_agent_pos!(agent, model)
  end

  return CellPopulation(model)
end

step(model::CellPopulation, stimulus::T) where T<:Stimulus = step(model, [stimulus])
step(model::CellPopulation, stimuli::Vector{T}=Vector{Stimulus}()) where T<:Stimulus = step!(model.model, (a, _) -> cellStepper(a, model, stimuli), modelStepper)

function cellStepper(agent::CellAgent, model::CellPopulation, stimuli::Vector{T}) where T<:Stimulus
  Δt = model.model.properties[:Δt]
  time = modelTime(model)
  doStep(agent, time, Δt, model, stimuli)
end


modelStepper(model::AgentBasedModel) = model.properties[:step_cnt] += 1

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
    die(cell)
    kill_agent!(agent, model.model)
  end

  for stimulus in stimuli
    stimulate(cell, stimulus, time)
  end

  if willDivide
    new_cell = divide(cell, time)
    if new_cell ≠ nothing
      new_agent = CellAgent(model.maxid[]+1, new_cell)
      add_agent_pos!(new_agent, model.model)
    end
  end
end

end
