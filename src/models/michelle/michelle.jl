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

include("MiCell.jl")

struct WildType <: CellType end
struct KO <: CellType end

#------------------ Cell factory --------------------
function cellFactory(parameters::Parameters, birth::Float64=0.0, cellType::T=GenericCell()) where T <: CellType
  michelleCell = Cell(birth, cellType)

  if T == WildType
    threshold = parameters.threshold
  else
    threshold = 0.0 # knockout nevers dies!
  end

  weights = copy(initialWeights)
  weights["BCLxL"] = parameters.bclxlWeight
  death = ThresholdDeath(threshold, weights, () -> TimeCourseParms(parameters.gstd))
  addTimer(michelleCell, death)

  division = DivisionTimer(λ_firstDivision, λ_divisionDestiny)
  addTimer(michelleCell, division)

  return michelleCell
end
#----------------------------------------------------

function run(model::CellPopulation, runDuration::Float64)
  Δt = modelTimeStep(model)
  time = 0:Δt:runDuration
  count = zeros(Int, length(time))
  cohort = zeros(Float64, length(time))

  proteinSampleTimes = Set([72.0, 96.0, 120.0, 140.0, 160.0, 180.0, 200.0])
  proteinLevels = DataFrame(time=Float64[], protein=String[], level=Float64[], genotype=String[])
  deathTimes = Float64[]
  sizehint!(deathTimes, length(model)*10)
  deathCounter(::Cell, time::Float64) = push!(deathTimes, time)
  model.deathCallback = deathCounter
  genotype = string(cellType(first(values(model.cells))))

  for (i, tm) in enumerate(time)
    step(model)
    count[i]  = cellCount(model)
    cohort[i] = cohortCount(model)

    if tm in proteinSampleTimes
      
      for protein in proteins
        for cell in model.cells
          level = proteinLevel(cell, protein, tm)
          push!(proteinLevels, (tm, protein, level, genotype))
        end
      end

      for cell in model.cells
        level = ensemble(cell, tm)
        push!(proteinLevels, (tm, "ensemble", level, genotype))
      end
    end
  end

  counts = DataFrame(time=time, count=count, cohort=cohort, genotype=genotype)
  result = Result(counts, proteinLevels, deathTimes)

  return result
end


thresholds    = 1.0:1.0:4.0
gstds         = 0.1:0.2:0.4
bclxllWeights = [0.1, 1, 10]
parameters = ConcreteParameters[]
for threshold in thresholds
  for gstd in gstds
    for weight in bclxllWeights
      p = ConcreteParameters(threshold, gstd, weight)
     push!(parameters, p)
    end
  end
end

results = Dict{ConcreteParameters, Result}()
nCells = 7477

struct Job
  parameter::Parameters
  result::Result
  lk::ReentrantLock
  function Job(parameter::Parameters)
    new(parameter, Result(), ReentrantLock())
  end
end
struct JobRef
  job::Job
  cellType::CellType
end
lk = ReentrantLock()
results = Dict{Parameters, Result}()

function runJob(jobRef::JobRef, results::Dict{Parameters, Result})
  cellType = jobRef.cellType
  job = jobRef.job
  parameter = job.parameter
  result = job.result

  model = createPopulation(nCells, (birth) -> cellFactory(parameter, birth, cellType))
  r = run(model, 200.0);
  lock(job.lk) do
    append!(result, r)
  end

  lock(lk) do 
    results[parameter] = result
  end
end

@info("Away we go!")

jobs = Vector{JobRef}()
for parameter in parameters
  job = Job(parameter)
  for cellType in [WildType(), KO()]
    jobRef = JobRef(job, cellType)
    push!(jobs, jobRef)
  end
end
  

result = Result()
@threads for job in jobs
  runJob(job, results)
  @info("$(job.job.parameter) $(job.cellType) done")
end

@info("Runs finished!")

open("results2.dat", "w") do io
  serialize(io, results)
end

@info("Data saved!")

# include("plotting.jl")
# @info("Data plotted!")

@info("You are awesome!")

