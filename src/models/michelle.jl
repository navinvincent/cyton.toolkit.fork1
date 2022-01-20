import Cyton: 
  LogNormalParms,
  DeathAccumulator,
  CellModule,
  DifferentationAccumulator,
  DivisionAccumulator,
  DeathAccumulator,
  SimpleDecay,
  shouldDie,
  shouldDifferentiate,
  shouldDivide,
  Cell,
  draw,
  addModule

proteins = ["BCL2", "BCLxL", "MCL1", "BIM"]

"A protein half life of around 24 hours"
d_λ = LogNormalParms(3.2, 0.4)
"An amount close to 1"
d_amount = LogNormalParms(0.3, 0.4)

struct ThresholdDeath <: DeathAccumulator
  "The threshold below which the cell will die"
  threshold::Float64
  "The protein level weights"
  weights::AbstractDict{AbstractString, Float64}
  "References to the current protein level"
  proteinLevels::AbstractDict{AbstractString, CellModule}

  function ThresholdDeath(proteinLevels::Dict{AbstractString, CellModule}) 
    return new(10.0,
    Dict{AbstractString, Float64}([
      "BCL2"=>5,
      "BCLxL"=>5,
      "MCL1"=>5,
      "BIM"=>-5,
    ]),
    proteinLevels
    )
  end
end

function shouldDie(deathAccumulator::ThresholdDeath, time::Float64)
  sum = 0.0
  for we in deathAccumulator.weights
    p = we.first
    w = we.second
    sum += w * deathAccumulator.proteinLevels[p].amount
  end

  return sum < deathAccumulator.threshold
end

"""
There are 3 differentention states
* Before first division
* Dividing
* Division destiny
"""
struct PreDivisionAccumulator <: DifferentationAccumulator
  start_diving::Float64
end
cellType(differentationAccumulator::PreDivisionAccumulator) = "preparing to divide"

struct PostDivisionAccumulator <: DifferentationAccumulator end
cellType(differentationAccumulator::PostDivisionAccumulator) = "division destiny"

struct DividingDivisionAccumulator <: DifferentationAccumulator
  stop_dividing::Float64
end
cellType(differentationAccumulator::DividingDivisionAccumulator) = "dividing"

"""
There are 2 division states
* Before first division and after division destiny. Cells don't divide in this state
* The period when cells divide
"""
struct NonDividingDivisionAccumulator <: DivisionAccumulator end
shouldDivide(divisionAccucumulator::NonDividingDivisionAccumulator, time::Float64) = false

struct SimpleDivisionAccumulator <: DivisionAccumulator
  next_division::Float64
end
shouldDivide(divisionAccucumulator::SimpleDivisionAccumulator, time::Float64) = time >= divisionAccucumulator.next_division

function cellFactory(birth::Float64=0.0) 
  michelleCell = Cell(birth)
  michelleCell.divisionAccumulator = nothing
  michelleCell.differentiationAccumulator = nothing

  for p in proteins
    addModule(michelleCell, p, SimpleDecay(draw(d_amount), draw(d_λ)))
  end

  michelleCell.deathAccumulator = ThresholdDeath(michelleCell.modules)

  return michelleCell
end