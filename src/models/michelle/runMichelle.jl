using Distributed, ClusterManagers

function freeWorkers()
  for i in workers() # Tear down any remaining workers
    if i > 1
      rmprocs(i)
    end
  end  
end

freeWorkers()

@info("Acquiring remote resources")
addprocs(SlurmManager(20), partition="regular", t="01:10:00", mem_per_cpu="2G", unbuffered="")

@everywhere include("michelle.jl")
@everywhere include("plotting.jl")

@everywhere struct SimJobResult <: JobResult
  parameterKey::ParameterKey
  result::Result
end

@everywhere function runJob(parms::ConcreteParameters)
  r = runModel(parms)
  k = parameterKey(parms)
  @info("runModel: $(parms)")
  SimJobResult(k, r)
end

@everywhere struct PlotJobResult <: JobResult end
@everywhere function runJob(parms::PlotParameters)
  bigPlot(parms.actualParameters, parms.result)
  @info("bigPlot: $(parms)")
  return PlotJobResult()
end

nJobs = length(parameters)

# Channels for sending jobs and returning results
jobChannel = RemoteChannel(() -> Channel{Parameters}(nJobs))
resultChannel = RemoteChannel(() -> Channel{JobResult}(nJobs))

@everywhere function doRunJob(jobChannel, resultChannel)
  @info("runJob started")
  while true
    parms = take!(jobChannel)
    result = runJob(parms)
    put!(resultChannel, result)
  end
end

# Push jobs into the job channel
@info("Load parameters into the queue")
for parameter in parameters
  put!(jobChannel, parameter)
end

@info("Launching the remote workers")
# remote_do(doRunJob, 1, jobChannel, resultChannel)
for w in workers()
  remote_do(doRunJob, w, jobChannel, resultChannel)
end

# Retrieve results and store in a dictionary
results = Dict{ParameterKey, Result}()
@info("Waiting for simulation results")
for _ in 1:nJobs
  jobResult = take!(resultChannel)
  local p = jobResult.parameterKey
  local r = jobResult.result
  if haskey(results, p)
    append!(results[p], r)
  else
    results[p] = r
  end
end

@info("Saving results")
open("results.dat", "w") do io
  serialize(io, results)
end

@info("Load plot jobs into queue")
for (parms, result) in results
  pp = PlotParameters(parms, result)
  put!(jobChannel, pp)
end

@info("Waiting for plot jobs to finish")
for _ in 1:length(results)
  take!(resultChannel)
end


@info("Freeing resources")
freeWorkers()

@info("You are awesome!")

