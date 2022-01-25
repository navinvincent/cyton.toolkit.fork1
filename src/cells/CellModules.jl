abstract type CellModule end
abstract type Accumulator <: CellModule end
abstract type DeathAccumulator <: Accumulator end
abstract type DivisionAccumulator <: Accumulator end
abstract type FateAccumulator <: Accumulator end

"Step the cell module forward by one time increment"
step(cellModule::CellModule, time::Float64, dt::Float64) = error("Not implemented")
"If not module, do nothing"
step(_::Nothing, time::Float64, Δt::Float64) = nothing

"This is called every time step. If it return true the cell is removed from the simmulation"
shouldDie(deathAccumulator::DeathAccumulator, time::Float64) = error("Not implemented")
shouldDie(_::Nothing, _::Float64) = false

"Call every time step and potentiall returns a new cell"
shouldDivide(divisionAccumulator::DivisionAccumulator, time::Float64) = error("Not implemented")
shouldDivide(_::Nothing, _::Float64) = false

"Called every time step to determine of the cell should differentiate"
shouldChangeFate(fateAccumulator::FateAccumulator, time::Float64) = error("Not implemented")
shouldChangeFate(_::Nothing, _::Float64) = false

"Provides a string to indicate the cell type"
cellType(differentationAccumulator::FateAccumulator) = "-- unknown --"

"A simple exponential decay of protein level"
mutable struct SimpleDecay <: Accumulator
  "Current amount of this protein"
  amount::Float64
  "Protein decay time constant"
  λ::Float64
end

step(death::SimpleDecay, time::Float64, Δt::Float64) = death.amount = death.amount * exp(-Δt/death.λ)
