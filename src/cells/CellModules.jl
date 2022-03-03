abstract type FateTimer end

export shouldDie, shouldDivide, inherit

"Step the cell module forward by one time increment"
step(::FateTimer, time::Float64, Δt::Float64) = error("stepTimer method not implemented")

"This is called every time step. If it return true the cell is removed from the simmulation"
shouldDie(::FateTimer, time::Float64) = false

"Call every time step and potentiall returns a new cell"
shouldDivide(::FateTimer, time::Float64) = false

"Inheritence mechanism for fate timers"
inherit(timer::FateTimer, time::Float64) = error("inherit method not implemented")

"A simple exponential decay of protein level"
mutable struct PoissonTimer <: FateTimer
  "Current amount of this stuff"
  amount::Float64
  "Decay time constant"
  λ::Float64
end

step(timer::PoissonTimer, _::Float64, Δt::Float64) = timer.amount = timer.amount * exp(-Δt/timer.λ)
