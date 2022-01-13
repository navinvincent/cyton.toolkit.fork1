using Agents, DistributionParms, CellModules

mutable struct Cell <: AbstractAgent
  "These are required by the ABM framework"
  id::Int
  pos::NTuple{2, Int}

  "Cell specific state"
  birth::Float64
  modules::AbstractArray{CellModule}
  deathAccumulator::DeathAccumulator
  divisionAccumulator::DivisionAccumulator
  
  function Cell(id::Int, pos::NTuple{2, Int}, birth::Float64)
    return new(id, pos, birth, [])
  end
end

add_module(cell::Cell, cellModule::CellModule) = append!(cell.modules, cellModule)

function step_cell(cell::Cell, time::Float64, model::AgentBasedModel)
  if time > cell.λ + cell.birth
    cell.birth = time
    cell.λ = draw(cell.d)
    cell.deaths += 1
  end
  for cellModule in cell.modules
    step_module(cellModule, time)
  end
  if should_die(cell.deathAccumulator)
    die(cell)
    kill_agent(cell, model)
  end
  if should_divide(cell.divisionAccumulator)
    new_cell = divide(cell)
    add_agent_single!(new_cell, model)
  end
end

age(cell::Cell, time::Float64) = time - cell.birth
remaining(cell::Cell, time::Float64) = cell.λ - age(cell, time)

die(cell::Cell) = nothing
divide(cell::Cell) = error("not implemented")
