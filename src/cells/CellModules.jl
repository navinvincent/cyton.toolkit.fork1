abstract type FateTimer end

export shouldDie, shouldDivide, inherit

"Step the cell module forward by one time increment"
step(timer::FateTimer, ::Time, ::Duration)::Union{CellEvent, Nothing} = error("step method not implemented for $(typeof(timer))")

"This is called every time step. If it return true the cell is removed from the simmulation"
shouldDie(::FateTimer, ::Time) = false

"Call every time step and potentiall returns a new cell"
shouldDivide(::FateTimer, ::Time) = false

"Inheritence mechanism for fate timers"
inherit(::FateTimer, ::Time) = error("inherit method not implemented for $(typeof(timer))")

"A simple exponential decay of protein level"
mutable struct PoissonTimer <: FateTimer
  "Current amount of this stuff"
  amount::Float64
  "Decay time constant"
  λ::Float64
end

step(timer::PoissonTimer, _::Time, Δt::Duration) = timer.amount = timer.amount * exp(-Δt/timer.λ)
