abstract type FateTimer end

"Step the cell module forward by one time increment"
step(_::FateTimer, time::Float64, Δt::Float64) = error("step method not implemented")

"This is called every time step. If it return true the cell is removed from the simmulation"
shouldDie(_::FateTimer, time::Float64) = false

"Call every time step and potentiall returns a new cell"
shouldDivide(fateTimer::FateTimer, time::Float64) = false

"Inheritence mechanism for fate timers"
inherit(fateTImer::FateTimer, time::Float64) = error("inherit method not implemented")

"A simple exponential decay of protein level"
mutable struct SimpleDecay <: FateTimer
  "Current amount of this protein"
  amount::Float64
  "Protein decay time constant"
  λ::Float64
end

step(death::SimpleDecay, _::Float64, Δt::Float64) = death.amount = death.amount * exp(-Δt/death.λ)
