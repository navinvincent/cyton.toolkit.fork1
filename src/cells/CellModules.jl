abstract type CellModule end
abstract type Accumulator <: CellModule end
abstract type DeathAccumulator <: Accumulator end
abstract type DivisionAccumulator <: Accumulator end

step_module(cellModule::CellModule, time::Float64) = nothing

should_die(deathAccumulator::DeathAccumulator) = false
should_divide(deathAccumulator::DeathAccumulator) = false
