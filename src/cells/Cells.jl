using Agents

abstract type AbstractCell end
abstract type Stimulus end

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
  divisonCount::Int64
  modules::AbstractDict{AbstractString, CellModule}
  deathAccumulator::Union{DeathAccumulator, Nothing}
  differentiationAccumulator::Union{DifferentationAccumulator, Nothing}
  divisionAccumulator::Union{DivisionAccumulator, Nothing}
  
  function Cell(birth::Float64)
    return new(birth, 0, Dict{AbstractString, CellModule}(), nothing, nothing, nothing)
  end

  function Cell(cell::Cell)
    return new(cell.birth, cell.divisonCount, Dict{AbstractString, CellModule}(), nothing, nothing, nothing)
  end

end

stimulate(cell::Cell, stimulus::Stimulus) = nothing
addModule(cell::Cell, name::AbstractString, cellModule::CellModule) = cell.modules[name] = cellModule

function step(agent::CellAgent, time::Float64, Δt::Float64, model::AgentBasedModel)
  cell = agent.cell

  for cellModule in cell.modules
    step(cellModule.second, time, Δt)
  end

  step(cell.deathAccumulator, time, Δt)
  step(cell.differentiationAccumulator, time, Δt)
  step(cell.divisionAccumulator, time, Δt)

  if shouldDie(cell.deathAccumulator, time)
    die(cell)
    kill_agent!(agent, model)
  end

  if shouldDivide(cell.divisionAccumulator, time)
    new_cell = divide(cell, time)
    new_agent = CellAgent(length(model.agents)+1, new_cell)
    add_agent_single!(new_agent, model)
  end

  if shouldDifferentiate(cell.differentiationAccumulator, time)
    differentiate(cell, time)
  end
end

age(cell::Cell, time::Float64) = time - cell.birth

die(cell::AbstractCell) = nothing
divide(cell::AbstractCell, time::Float64) = error("Not implemented")
differentiate(cell::AbstractCell, time::Float64) = error("Not implemented")
