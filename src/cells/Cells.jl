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
  generation::Int64
  timers::Vector{FateTimer}
  
  function Cell(birth::Float64)
    return new(birth, 0, [])
  end

  function Cell(birth::Float64, divisionCount::Int64)
    return new(birth, divisionCount, [])
  end

end

stimulate(_::Cell, _::Stimulus) = nothing
addTimer(cell::Cell, timer::FateTimer) = push!(cell.timers, timer)

function step(agent::CellAgent, time::Float64, Δt::Float64, model::AgentBasedModel)
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
    new_agent = CellAgent(model.maxid[]+1, new_cell)
    add_agent_pos!(new_agent, model)
  end

end

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

