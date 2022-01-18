using Agents

mutable struct CellAgent <: AbstractAgent
  "These are required by the ABM framework"
  id::Int
  pos::NTuple{2, Int}
  cell::Cell
end

abstract type AbstractCell end
struct Cell <: AbstractCell
  "Cell specific state"
  birth::Float64
  modules::AbstractDict{AbstractString, CellModule}
  deathAccumulator::DeathAccumulator
  diffentiationAccumulator::DiffentiationAccumulator
  divisionAccumulator::DivisionAccumulator
  
  function Cell(birth::Float64)
    return new(birth, Dict{String, CellModule}(), nothing, nothing, nothing)
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
    kill_agent(cell, model)
  end

  if shouldDivide(cell.divisionAccumulator, time)
    new_cell = divide(cell, time)
    add_agent_single!(new_cell, model)
  end

  if shouldDifferentiate(cell.diffentationAccumulator, time)
    differentiate(cell, time)
  end
end

age(cell::Cell, time::Float64) = time - cell.birth

die(cell::Cell) = nothing
divide(cell::Cell, time::Float64) = error("not implemented")
differentiate(cell::Cell, time::Float64) = error("not implemented")
