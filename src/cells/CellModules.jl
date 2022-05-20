abstract type FateTimer end

export shouldDie, shouldDivide, inherit

"Step the cell module forward by one time increment"
step(timer::FateTimer, time::Time, Δt::Duration)::Union{CellEvent, Nothing} = error("step method not implemented for $(typeof(timer))")

"This is called every time step. If it return true the cell is removed from the simmulation"
shouldDie(::FateTimer, time::Time) = false

"Call every time step and potentiall returns a new cell"
shouldDivide(::FateTimer, time::Time) = false

"Inheritence mechanism for fate timers"
inherit(timer::FateTimer, time::Time) = error("inherit method not implemented")

"A simple exponential decay of protein level"
mutable struct PoissonTimer <: FateTimer
  "Current amount of this stuff"
  amount::Float64
  "Decay time constant"
  λ::Float64
end

step(timer::PoissonTimer, _::Float64, Δt::Duration) = timer.amount = timer.amount * exp(-Δt/timer.λ)
