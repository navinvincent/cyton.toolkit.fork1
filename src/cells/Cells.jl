using Agents

abstract type AbstractCell end

mutable struct CellAgent <: AbstractAgent
  "These are required by the ABM framework"
  id::Int
  pos::NTuple{2, Int}
  cell::AbstractCell
  function CellAgent(cell::AbstractCell)
    return new(0, (0, 0), cell)
  end
end

mutable struct Cell <: AbstractCell
  birth::Float64
  modules::AbstractDict{AbstractString, CellModule}
  deathAccumulator::Union{DeathAccumulator, Nothing}
  differentiationAccumulator::Union{DifferentationAccumulator, Nothing}
  divisionAccumulator::Union{DivisionAccumulator, Nothing}
  
  function Cell(birth::Float64)
    return new(birth, Dict{AbstractString, CellModule}(), nothing, nothing, nothing)
  end
end

addModule(cell::Cell, name::AbstractString, cellModule::CellModule) = cell.modules[name] = cellModule

function step(agent::CellAgent, time::Float64, dt::Float64, model::AgentBasedModel)
  cell = agent.cell

  for cellModule in cell.modules
    step(cellModule.second, time, dt)
  end

  if shouldDie(cell.deathAccumulator, time)
    die(cell)
    kill_agent!(agent, model)
  end

  if shouldDivide(cell.divisionAccumulator, time)
    new_cell = divide(cell, time)
    new_agent = CellAgent(new_cell)
    add_agent_single!(new_agent, model)
  end

  if shouldDifferentiate(cell.differentiationAccumulator, time)
    differentiate(cell, time)
  end
end

age(cell::Cell, time::Float64) = time - cell.birth

die(cell::Cell) = nothing
divide(cell::Cell, time::Float64) = error("Not implemented")
differentiate(cell::Cell, time::Float64) = error("Not implemented")
