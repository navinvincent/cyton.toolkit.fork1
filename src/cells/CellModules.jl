abstract type CellModule end
abstract type Accumulator <: CellModule end
abstract type DeathAccumulator <: Accumulator end
abstract type DivisionAccumulator <: Accumulator end
abstract type DiffentationAccumulator <: Accumulator end

"Step the cell module forward by one time increment"
step(cellModule::CellModule, time::Float64, dt::Float64) = error("Not implemented")

"This is called every time step. If it return true the cell is removed from the simmulation"
shouldDie(deathAccumulator::DeathAccumulator, time::Float64) = error("Not implemented")

"Call every time step and potentiall returns a new cell"
shouldDivide(divisionAccumulator::DivisionAccumulator, time::Float64) = error("Not implemented")

"Called every time step to determine of the cell should differentiate"
shouldDifferentiate(diffentationAccumulator::DiffentationAccumulator, time::Float64) = error("Not implemented")

"Provides a string to indicate the cell type"
cellType(diffentationAccumulator::DiffentationAccumulator) = "-- unknown --"

"A simple exponential decay of protein level"
mutable struct SimpleDecay <: CellModule
  "Current amount of this protein"
  amount::Float64
  "Protein half life"
  λ::Float64
end

step(death::SimpleDecay, time::Float64, dt::Float64) = death.amount = death.amount = dt*death.λ
