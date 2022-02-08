module Cyton

using Agents

include("probability/DistributionParms.jl")
include("cells/CellModules.jl")
include("cells/Cells.jl")
include("utils/void_space.jl")


model_time(model::AgentBasedModel) = model.properties[:step_cnt] * model.properties[:Δt]

function createModel(Ncells::Int, cellFactory)
  space = VoidSpace()
  scheduler = Schedulers.fastest
  properties = Dict(:step_cnt => 0, :Δt => 0.1)
  model = AgentBasedModel(CellAgent, space; properties, scheduler)

  for id in 1:Ncells
    cell = cellFactory()
    agent = CellAgent(id, cell)
    add_agent_pos!(agent, model)
  end

  return model
end

step(model::AgentBasedModel) = step!(model, cellStepper, modelStepper)

function cellStepper(agent::CellAgent, model::AgentBasedModel)
  Δt = model.properties[:Δt]
  time = model_time(model)
  doStep(agent, time, Δt, model)
end

modelStepper(model::AgentBasedModel) = model.properties[:step_cnt] += 1

function doStep(agent::CellAgent, time::Float64, Δt::Float64, model::AgentBasedModel)
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
    kill_agent!(agent, model)
  end

  if willDivide
    new_cell = divide(cell, time)
    if new_cell ≠ nothing
      new_agent = CellAgent(model.maxid[]+1, new_cell)
      add_agent_pos!(new_agent, model)
    end
  end

end

end
