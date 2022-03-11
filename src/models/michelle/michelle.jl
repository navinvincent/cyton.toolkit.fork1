"""
A cell model based on a mix of Michelle's experiments and the Cyton2 paper.
- Time to first division and destiny are based on the Cyton2 data
- Death is based on a model of protein levels over time. An ensemble value
  is calculated and if the ensemble is less than a threshold the cell dies.
  The ensemble is a weighted sum of the individual protein levels. 

"""

using Cyton
import Cyton: shouldDie, shouldDivide, inherit, step

using DataFrames, Colors, Cairo

using Base.Threads, Serialization
import Base.show

# Gadfly defaults
Gadfly.set_default_plot_size(20cm, 20cm)
Gadfly.push_theme(Theme(background_color="white"))

include("MiCell.jl")

#------------------ Cell factory --------------------
function cellFactory(parameters::Parameters, birth::Float64=0.0)
  michelleCell = Cell(birth)

  death = ThresholdDeath(parameters.threshold, initialWeights, TimeCourseParms(parameters.gstd))
  addTimer(michelleCell, death)
  division = DivisionTimer(λ_firstDivision, λ_divisionDestiny)
  addTimer(michelleCell, division)

  return michelleCell
end
#----------------------------------------------------

struct Result
  counts::DataFrame
  proteinLevels::DataFrame
  deathTimes::Vector{Float64}
end

function run(model::CellPopulation, runDuration::Float64)
  Δt = modelTimeStep(model)
  time = 0:Δt:runDuration
  count = zeros(Int, length(time))

  proteinSampleTimes = Set([72.0, 96.0, 120.0, 144.0, 168.0])
  proteinLevels = DataFrame(time=Float64[], name=String[], level=Float64[])
  deathTimes = Float64[]
  deathCounter(::Cell, time::Float64) = push!(deathTimes, time)
  model.deathCallback = deathCounter

  for (i, tm) in enumerate(time)
    step(model)
    count[i] = cellCount(model)

    if tm in proteinSampleTimes
      
      for protein in proteins
        for cell in model.cells
          level = proteinLevel(cell, protein, tm)
          push!(proteinLevels, (tm, protein, level))
        end
      end

      for cell in model.cells
        level = ensemble(cell, tm)
        push!(proteinLevels, (tm, "ensemble", level))
      end
    end
  end

  counts = DataFrame(time=time, count=count)
  result = Result(counts, proteinLevels, deathTimes)

  return result
end

struct ConcreteParameters <: Parameters
  threshold::Float64
  gstd::Float64
end
show(io::IO, ::Parameters) = print(io, "no desciption")
show(io::IO, parms::ConcreteParameters) = print(io, "threshold=$(parms.threshold) gstd=$(parms.gstd)")

lk = ReentrantLock()
results = Dict{Parameters, Result}()

thresholds = 1.0:2.0:7.0
gstds      = 0.1:0.2:0.5
parameters = ConcreteParameters[]
for threshold in thresholds
  for gstd in gstds
    p = ConcreteParameters(threshold, gstd)
    push!(parameters, p)
  end
end

@info("Away we go!")
results = Dict{ConcreteParameters, Result}()
@threads for parameter in parameters
  model = createPopulation(7477, (birth) -> cellFactory(parameter, birth))
  result = run(model, 200.0);
  lock(lk) do 
    results[parameter] = result
  end
  @info("$(parameter) done")
end

@info("Runs finished!")

open("results.dat", "w") do io
  serialize(io, results)
end

@info("Data saved!")

include("plotting.jl")

@info("Data plotted!")

@info("You are awesome!")

