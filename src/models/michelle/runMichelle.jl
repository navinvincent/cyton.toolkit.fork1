using Distributed, ClusterManagers

@info("Acquiring remote resources")
addprocs(SlurmManager(2), partition="regular", t="01:10:00", mem_per_cpu="2G")

@everywhere include("michelle.jl")
@everywhere include("plotting.jl")

function f()
  flush(stdout)
  flush(stderr)
end

struct SimJobResult <: JobResult
  parameterKey::ParameterKey
  result::Result
end

function runJob(parms::ConcreteParameters)
  r = runModel(parms)
  k = parameterKey(parms)
  @info("runModel: $(parms)")
  f()
  SimJobResult(k, r)
end

struct PlotJobResult <: JobResult end
function runJob(parms::PlotParameters)
  bigPlot(parms.actualParameters, parms.result)
  @info("bigPlot: $(parms)")
  f()
  return PlotJobResult()
end


# Channels for sending jobs and returning results
jobChannel = RemoteChannel(() -> Channel{Parameters}(50))
resultChannel = RemoteChannel(() -> Channel{JobResult}(50))

@everywhere function doRunJob(jobChannel, resultChannel)
  @info("runJob started")
  f()
  while true
    parms = take!(jobChannel)
    t = Task(() -> runJob(parms))
    schedule(t)
    errormonitor(t)
    result = fetch(t)
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
for _ in 1:length(parameters)
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
for i in workers()
  rmprocs(i)
end  

@info("You are awesome!")

