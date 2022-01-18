using Cyton

proteins = ["BCL2", "BCLxL", "MCL1", "BIM"]

"A protein half life of around 24 hours"
d_λ = LogNormalParms(3.2, 0.4)
"An amount close to 1"
d_amount = LogNormalParms(0.3, 0.4)

struct ThresholdDeath <: DeathAccumulator
  "The threshold below which the cell will die"
  threshold::Float64
  "The protein level weights"
  weights::Dict{AbstractString, Float64}
  "References to the current protein leveld"
  proteins::Dict{AbstractString, SimpleDecay}

  function ThresholdDeath(proteins::Dict{AbstractString, SimpleDecay}) 
    return new(10.0,
    weights = Dict{AbstractString, Float64}([
      "BCL2"=>5,
      "BCLxL"=>5,
      "MCL1"=>5,
      "BIM"=>-5,
    ]),
    proteins
    )
  end
end

function shouldDie(deathAccumulator::ThresholdDeath, time::Float64)
  sum = 0.0
  for we in deathAccumulator.weights
    p = we.first
    w = w.second
    sum += w * proteins[p]
  end

  return sum < deathAccumulator.threshold
end


"""
There are 3 differentention states
* Pre first division
* Dividing
* Division destiny
"""
struct PreDivisionAccumulator <: DiffentationAccumulator 
  start_diving::Float64
end
cellType(diffentationAccumulator::PreDivisionAccumulator) = "preparing to divide"

struct PostDivisionAccumulator <: DiffentationAccumulator end
cellType(diffentationAccumulator::PostDivisionAccumulator) = "division destiny"

struct DividingDivisionAccumulator <: DiffentationAccumulator
  stop_dividing::Float64
end
cellType(diffentationAccumulator::DividingDivisionAccumulator) = "dividing"

"""
There are 2 division states
* Before first division and after division destiny. Cells don't divide in this state
* The perion when cells divide
"""
struct NonDividingDivisionAccumulator <: DivisionAccumulator end
shouldDivide(divisionAccucumulator::NonDividingDivisionAccumulator, time::Float64) = false

struct SimpleDivisionAccumulator <: DivisionAccumulator
  next_division::Float64
end
shouldDivide(divisionAccucumulator::SimpleDivisionAccumulator, time::Float64) = time >= divisionAccucumulator.next_division

function createMichelleCell(birth::Float64=0.0) 
  michelleCell = Cell(birth)

  for p in proteins
    addModule(p, SimpleDecay(draw(d_amount), draw(d_λ)))
  end

  michelleCell.deathAccumulator = ThresholdDeath(michelleCell.modules)

  return michelleCell
end